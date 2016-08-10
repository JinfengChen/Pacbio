#!/usr/bin/perl
=header
Map Pacbio reads bax.h5 reads to assembly. Generate cmp.h5 and bam file. cmp.h5 can be used to polish assembly using quiver.  
USE fullpath of files.
--ref:   reference sequence
--input: fofn file which list all SMRT cell each include *.p0.1/2/3.bax.h5
--tool:  mapping tools: pbalign
--split: how many sequence in each sub fasta file
-project: project name that used for result file 

perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/yeast_ass.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.yeast.fofn --output yeast_ass_reads --project yeast_ass_reads --verbose > log 2>&1 &
=cut

use Getopt::Long;
use File::Basename;
use File::Spec;
use FindBin qw ($Bin);


my %opt;
GetOptions(\%opt,"ref:s","input:s","tool:s","cpu:s","bam","verbose","output:s", "split:s", "step:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref fullpath/all.con -1 fullpath/1.fq -2 fullpath/2.fq -min 0 -max 500 -cpu 12 -tool bwa\n";
   exit();
}


$opt{tool} ||= "pbalign";
$opt{cpu} ||=12;
$opt{split} ||= 2000;

unless ($opt{project}){
    print "Need to specify project name: current folder/out\n";
    exit(2)
}

$opt{output} = File::Spec->rel2abs("$opt{project}");
my $tempdir = File::Spec->rel2abs("./tempdir");
`mkdir $opt{output}` unless (-d $opt{output});
`mkdir "$opt{output}\_cmp_split"` unless (-d "$opt{output}\_cmp_split");

my $bwa="/opt/tyler/bin/";
my $blasr="/opt/blasr/453c25ab/bin//blasr";
my $soap="/usr/local/bin";
my $ssaha="/home/jfchen/software/ssaha2_v2.5.5_x86_64/";
my $maq="/opt/tyler/bin/maq";
my $pbalign="~/BigData/00.RD/Assembly/Pacbio/install/pythonlib/bin/pbalign";
my $bax2bam="/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/pitchfork/workspace/blasr/utils/bax2bam/bin/bax2bam";
my $python_bin="/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/pythonlib/bin";
my $fqsplit="/rhome/cjinfeng/software/bin/fastq_split.pl";
my $fasplit="$Bin/fastaDeal.pl";
my $SAMtool="/opt/tyler/bin/samtools";
my $rmdup="/opt/picard/1.81/MarkDuplicates.jar";
my $loadchemistry='/opt/linux/centos/7.x/x86_64/pkgs/python/2.7.5/bin/loadChemistry.py';
my $loadpulses='/opt/blasr/453c25ab/bin/loadPulses';

if (exists $opt{input}){
   if ($opt{tool}=~/pbalign/){
      print "Run long read mapping by blasr!\n";
      ### split into small files
      my @split;
      my @h5files = readfofn($opt{input}, $opt{output});

      ### map small files
      #m130605_000141_42207_c100515142550000001823076608221372_s1_p0.1.bax.h5
      my @map;
      my @fqs = @h5files;
      my @cmph5;
      my @samdirs;
      for(my $i=0; $i<@fqs; $i+=1){
           print "$i\t$fqs[$i]\n";
           #my $prefix=substr(basename($fqs[$i]),0,3);
           my $temp=basename($fqs[$i]);
           my $prefix=$1 if ($temp=~/^(.*)\.fofn/);
           my $temp_sam="$opt{output}/$prefix.tempdir";
           #new pipeline need pbalign to have bam in and bam out: https://github.com/PacificBiosciences/pbalign/wiki/Mapping-RS-II-reads-(in-movie.bax.h5)-and-calling-consensus-using-pbalign-BLASR-while--sam-option-is-deprecated
           my @bax_files = read_baxh5($fqs[$i]);
           push @map, "$bax2bam $bax_files[0] $bax_files[1] $bax_files[2] -o $opt{output}/$prefix --subread";
           push @map, "$pbalign --algorithm blasr --tmpDir $temp_sam --nproc 16 $opt{output}/$prefix.subreads.bam $opt{ref} $opt{output}/$prefix.bam";
           push @cmph5, "$opt{output}/$prefix.bam";
           push @samdirs, $temp_sam;
           #push @map, "$blasr $fqs[$i] $opt{ref} -sam -bestn 2 -nproc 1 > $opt{output}/$prefix.sam";
      }
      my $cmd1=join("\n",@map);
      writefile("$opt{project}.map.sh","$cmd1\n");
      if ($opt{step}=~/1/){
          `perl qsub-pbs-env_bam.pl --maxjob 40 --lines 2 --interval 120 --resource nodes=1:ppn=16,walltime=40:00:00,mem=40g --convert no $opt{project}.map.sh`;
      }

      ### merge and clean tmp files
      my @merge;
      my $cmph5_files = join(" ", @cmph5);
      push @merge, "python $python_bin/cmph5tools.py merge --outFile $opt{output}.cmp.h5 $cmph5_files" unless (-e "$opt{output}.cmp.h5");
      #push @merge, "python $python_bin/cmph5tools.py sort --deep $opt{output}.cmp.h5 --outFile $opt{output}.cmp.dsort.h5 --tmpDir $tempdir" unless (-e "$opt{output}.cmp.dsort.h5");
      #push @merge, "python $python_bin/cmph5tools.py select --groupBy=Reference --outDir $opt{output}\_cmp_split $opt{output}.cmp.h5";
      if ($opt{bam}){
          my @samfiles = getsamfiles(\@samdirs);
          unless (-e "$opt{output}.sam"){
              push (@merge, "grep \"^@\" $samfiles[0] > $opt{output}.sam");
              for (my $i=0; $i<@samfiles;$i++){
                  push (@merge, "cat $samfiles[$i] | grep -v \"^@\" >> $opt{output}.sam");
              }
          }
          push (@merge, "$SAMtool view -bS -o $opt{output}.raw.bam $opt{output}.sam > $opt{output}.convert.log 2> $opt{output}.convert.log2") unless (-e "$opt{output}.raw.bam");
          push (@merge, "$SAMtool sort -m 1000000000 $opt{output}.raw.bam $opt{output} > $opt{output}.sort.log 2> $opt{output}.sort.log2") unless (-e "$opt{output}.bam");
      }
      #push (@merge, "java -Xmx5G -jar $rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$opt{output}.sort.bam OUTPUT=$opt{output}.bam METRICS_FILE=$opt{output}.dupli > $opt{output}.rmdup.log 2> $opt{output}.rmdup.log2") unless (-e "$opt{output}.bam");
      #push (@merge, "$SAMtool index $opt{output}.bam");
      my $cmd2=join("\n",@merge);
      writefile("$opt{project}.merge.sh","$cmd2\n");
      if ($opt{step}=~/2/){
          `perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs-env.pl --lines 5000 --interval 120  --resource nodes=1:ppn=1,walltime=100:00:00,mem=10G --convert no $opt{project}.merge.sh`;
      }

      my $ref_len = getfastalen($opt{ref});
      my @ref_id  = sort {$a <=> $b} keys %$ref_len;
      ###split and sort cmp.h5 files
      my @sort_split;
      for (my $i=0; $i<@ref_id; $i++){
          push @sort_split, "python $python_bin/cmph5tools.py select --where \"Reference == '$ref_id[$i]'\" --outFile $opt{output}\_cmp_split/$ref_id[$i].cmp.h5 $opt{output}.cmp.h5";
          push @sort_split, "python $python_bin/cmph5tools.py sort --deep --inPlace $opt{output}\_cmp_split/$ref_id[$i].cmp.h5 --tmpDir $tempdir";
          #cmph5tools in path of qsub, no need to do load* here
          #push @sort_split, "$loadpulses $opt{input} $opt{output}\_cmp_split/$ref_id[$i].cmp.h5 -metrics DeletionQV,DeletionTag,InsertionQV,MergeQV,SubstitutionQV";
          #push @sort_split, "$loadchemistry $opt{input} $opt{output}\_cmp_split/$ref_id[$i].cmp.h5";
          push @sort_split, "/usr/bin/h5repack -f GZIP=1 $opt{output}\_cmp_split/$ref_id[$i].cmp.h5 $opt{output}\_cmp_split/$ref_id[$i].cmp.h5.TMP";
          push @sort_split, "mv $opt{output}\_cmp_split/$ref_id[$i].cmp.h5.TMP $opt{output}\_cmp_split/$ref_id[$i].cmp.h5";
      }

      #my @cmph5_split_files = glob("$opt{output}\_cmp_split/*.cmp.h5");
      #my @sort_split;
      #for (my $i=0; $i<@cmph5_split_files; $i++){
      #    my $contig = $cmph5_split_files[$i] =~ /(.*)\.cmp\.h5/ ? $1 : 'NA';
      #    push @sort_split, "python $python_bin/cmph5tools.py sort --deep $contig.cmp.h5 --outFile $contig.sort.cmp.h5 --tmpDir $tempdir";
      #    #cmph5tools in path of qsub, no need to do load* here
      #    #push @sort_split, "$loadpulses $opt{input} $contig.sort.cmp.h5 -metrics DeletionQV,DeletionTag,InsertionQV,MergeQV,SubstitutionQV";
      #    #push @sort_split, "$loadchemistry $opt{input} $contig.sort.cmp.h5";
      #} 
      my $cmd3=join("\n", @sort_split);
      writefile("$opt{project}.sort_split.sh","$cmd3\n");
      if ($opt{step}=~/3/){
          #lines=60 need be dividable by 4 and 6 as we need 6 line for load and 4 line without load
          #for 3000 contig lines=60 will have 10 contig in each job and 300 jobs in total 
          `perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs-env.pl --maxjob 100 --lines 60 --interval 120  --resource nodes=1:ppn=1,walltime=100:00:00,mem=10G --convert no $opt{project}.sort_split.sh`;
      }

      ### clear tmp files
      my @clear;
      push @clear, "rm $opt{output}.sam $opt{output}.raw.bam";
      push @clear, "rm $opt{output}.sort.* $opt{output}.convert.*";
      push @clear, "rm $opt{output} $opt{output}.map.sh* $opt{output}.merge.sh* -R";
      push @clear, "rm $opt{output}.sort_split.sh* -R";
      push @clear, "rm $opt{output} -R";
      my $cmd4=join("\n",@clear);
      writefile("$opt{project}.clear.sh","$cmd4\n");
      unless ($opt{verbose}){
          `perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --lines 10 --interval 120 --convert no $opt{project}.clear.sh`;
      }
   }
}



sub writefile
{
my ($file,$line)=@_;
open WR, ">$file" or die "$!";
     print WR "$line";
close WR;
}

#in each sam directory get the sam file with smaller size, which is filtered
sub getsamfiles
{
my ($subdir) = @_;
my @bams;
my %hash;
for (my $i=0;$i<@$subdir;$i++){
    my @tempbam = glob("$subdir->[$i]/*.sam");
    for (my $j=0;$j<@tempbam;$j++){
        $hash{$tempbam[$j]} = -s $tempbam[$j];
    }
    my @sortbam = sort {$hash{$a} <=> $hash{$b}} keys %hash;
    print $sortbam[0], "\t", $hash{$sortbam[0]}, "\n";
    print $sortbam[1], "\t", $hash{$sortbam[1]}, "\n";
    push @bams, $sortbam[0];
}
return @bams;
}

#read inputfofn file, which is a list of multi SMRT cell file including three bax.h5 file of each SMRT cell
#this function split each SMRT into one input.fofn as input for pbalign. ruturn array with all the fofn files
#*.1.bax.h5
#*.2.bax.h5
#*.3.bax.h5
sub readfofn
{
my ($file, $prefix)=@_;
my @subfofn;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    if ($unit[0]=~/\/(.*?)\.\d+\.bax\.h5/){
    #if ($unit[0]=~/\/(.*?)\.bas\.h5/){
        push @{$hash{$1}}, $unit[0];
    }
    print "$unit[0]\n";
    #$hash{$unit[0]}=1;
}
close IN;

foreach my $pre (sort keys %hash){
    my $subfiles = join("\n", @{$hash{$pre}});
    my $subreads = basename($pre);
    print "$subreads\n$subfiles\n";
    open SUBFILE, ">$prefix.$subreads.fofn";
    print SUBFILE "$subfiles\n";
    close SUBFILE;
    push @subfofn, "$prefix.$subreads.fofn";
}
return @subfofn;
}


sub getfastalen
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    $hash{$head}= length $seq;
}
close IN;
$/="\n";
return \%hash;
}




sub read_baxh5
{
my ($file)=@_;
my @file;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    chomp $_;
    push @file, $_
}
close IN;
return @file;
}

