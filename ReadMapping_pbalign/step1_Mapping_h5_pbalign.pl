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
GetOptions(\%opt,"ref:s","input:s","tool:s","cpu:s","bam","verbose","output:s", "split:s", "project:s","help");

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


my $bwa="/opt/tyler/bin/";
my $blasr="/opt/blasr/453c25ab/bin//blasr";
my $soap="/usr/local/bin";
my $ssaha="/home/jfchen/software/ssaha2_v2.5.5_x86_64/";
my $maq="/opt/tyler/bin/maq";
my $pbalign="/opt/linux/centos/7.x/x86_64/pkgs/python/2.7.5/bin/pbalign";
my $python_bin="/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/pythonlib/bin";
my $fqsplit="/rhome/cjinfeng/software/bin/fastq_split.pl";
my $fasplit="$Bin/fastaDeal.pl";
my $SAMtool="/opt/tyler/bin/samtools";
my $rmdup="/opt/picard/1.81/MarkDuplicates.jar";

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
           push @map, "$pbalign --forQuiver --tmpDir $temp_sam --nproc 4 $fqs[$i] $opt{ref} $opt{output}/$prefix.cmp.h5";
           push @cmph5, "$opt{output}/$prefix.cmp.h5";
           push @samdirs, $temp_sam;
           #push @map, "$blasr $fqs[$i] $opt{ref} -sam -bestn 2 -nproc 1 > $opt{output}/$prefix.sam";
      }
      my $cmd1=join("\n",@map);
      writefile("$opt{project}.map.sh","$cmd1\n");
      `perl /rhome/cjinfeng/software/bin/qsub-pbs-env.pl --maxjob 40 --lines 1 --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=30g --convert no $opt{project}.map.sh`;

      ### merge and clean tmp files
      my @merge;
      my $cmph5_files = join(" ", @cmph5);
      push @merge, "python $python_bin/cmph5tools.py merge --outFile $opt{output}.cmp.h5 $cmph5_files" unless (-e "$opt{output}.cmp.h5");
      push @merge, "python $python_bin/cmph5tools.py sort --deep $opt{output}.cmp.h5 --outFile $opt{output}.cmp.dsort.h5 --tmpDir $tempdir" unless (-e "$opt{output}.cmp.dsort.h5");
      my @samfiles = getsamfiles(\@samdirs);
      unless (-e "$opt{output}.sam"){
          push (@merge, "grep \"^@\" $samfiles[0] > $opt{output}.sam");
          for (my $i=0; $i<@samfiles;$i++){
              push (@merge, "cat $samfiles[$i] | grep -v \"^@\" >> $opt{output}.sam");
          }
      }
      push (@merge, "$SAMtool view -bS -o $opt{output}.raw.bam $opt{output}.sam > $opt{output}.convert.log 2> $opt{output}.convert.log2") unless (-e "$opt{output}.raw.bam");
      push (@merge, "$SAMtool sort -m 1000000000 $opt{output}.raw.bam $opt{output} > $opt{output}.sort.log 2> $opt{output}.sort.log2") unless (-e "$opt{output}.bam");
      #push (@merge, "java -Xmx5G -jar $rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$opt{output}.sort.bam OUTPUT=$opt{output}.bam METRICS_FILE=$opt{output}.dupli > $opt{output}.rmdup.log 2> $opt{output}.rmdup.log2") unless (-e "$opt{output}.bam");
      #push (@merge, "$SAMtool index $opt{output}.bam");
      my $cmd2=join("\n",@merge);
      writefile("$opt{project}.merge.sh","$cmd2\n");
      `perl /rhome/cjinfeng/software/bin/qsub-pbs-env.pl --lines 5000 --interval 120  --resource walltime=100:00:00,mem=10G --convert no $opt{project}.merge.sh`;

      ### clear tmp files
      my @clear;
      push @clear, "rm $opt{output}.sam $opt{output}.raw.bam";
      push @clear, "rm $opt{output}.sort.* $opt{output}.convert.*";
      push @clear, "rm $opt{output} $opt{output}.map.sh* $opt{output}.merge.sh* -R";
      my $cmd3=join("\n",@clear);
      writefile("$opt{project}.clear.sh","$cmd3\n");
      unless ($opt{verbose}){
          `perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --lines 3 --interval 120 --convert no $opt{project}.clear.sh`;
      }
   }elsif($opt{tool}=~/soap/){ # soap
      print "Run pair-end mapping by soap!\n";
      unless (-e "$opt{ref}.index.sai"){
         `$soap/2bwt-builder $opt{ref} > $opt{project}.builder.log 2> $opt{project}.builder.log2`;
         `$SAMtool faidx $opt{ref}`; # generate $opt{ref}.fai used in samtools view
      }
      `$soap/soap -a $opt{1} -b $opt{2} -D $opt{ref}.index -o $opt{project}.soap.PE -2 $opt{project}.soap.SE -p $opt{cpu} -m $opt{min} -x $opt{max} > $opt{project}.soap.log 2> $opt{project}.soap.log2` if ($opt{max} < 2000);
      `$soap/soap -a $opt{1} -b $opt{2} -D $opt{ref}.index -o $opt{project}.soap.PE -2 $opt{project}.soap.SE -p $opt{cpu} -m $opt{min} -x $opt{max} -R > $opt{project}.soap.log 2> $opt{project}.soap.log2` if ($opt{max} >= 2000);
      #`cat $opt{project}.soap.PE $opt{project}.soap.SE > $opt{project}.soap`;
      if ($opt{bam}){
      print "Convert SOAP to SAM\n";
      `perl /opt/tyler/bin/soap2sam.pl $opt{project}.soap.SE > $opt{project}.soap.SE.sam`;
      `perl /opt/tyler/bin/soap2sam.pl -p $opt{project}.soap.PE > $opt{project}.soap.PE.sam`;
      print "Convert SAM to BAM, sort and merge\n";
      `$SAMtool view -bS -t $opt{ref}.fai -o $opt{project}.raw.SE.bam $opt{project}.soap.SE.sam > $opt{project}.convert.SE.log 2> $opt{project}.convert.SE.log2`;
      `$SAMtool sort $opt{project}.raw.SE.bam $opt{project}.SE.sort > $opt{project}.SE.sort.log 2> $opt{project}.SE.sort.log2`;
      `$SAMtool view -bS -t $opt{ref}.fai -o $opt{project}.raw.PE.bam $opt{project}.soap.PE.sam > $opt{project}.convert.PE.log 2> $opt{project}.convert.PE.log2`;
      `$SAMtool sort $opt{project}.raw.PE.bam $opt{project}.PE.sort > $opt{project}.PE.sort.log 2> $opt{project}.PE.sort.log2`;
      `$SAMtool merge -f $opt{project}.sort.bam $opt{project}.SE.sort.bam $opt{project}.PE.sort.bam`;  
      print "Remove duplicate!\n";
      `java -jar $rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$opt{project}.sort.bam OUTPUT=$opt{project}.bam METRICS_FILE=$opt{project}.dupli > $opt{project}.rmdup.log 2> $opt{project}.rmdup.log2`;
      }
      unless ($opt{verbose}){
          `rm $opt{project}.*.sam $opt{project}.raw.*.bam $opt{project}.*.sort.bam $opt{project}.sort.bam`;
          `rm $opt{project}.*.log* $opt{project}.dupli`;
          
      }
      print "Done!\n";
   }elsif($opt{tool}=~/maq/){ # maq
      ## build reference index
      unless (-e "$opt{ref}.bfa"){
         `$maq fasta2bfa $opt{ref} $opt{ref}.bfa`;
      }
      `$maq fastq2bfq $opt{1} $opt{1}.bfq` unless (-e "$opt{1}.bfq");
      `$maq fastq2bfq $opt{2} $opt{2}.bfq` unless (-e "$opt{2}.bfq");
      `$maq match -a $opt{max} $opt{project}.Maq.map $opt{ref}.bfa $opt{1}.bfq $opt{2}.bfq`;
      `$maq mapview $opt{project}.Maq.map > $opt{project}.Maq.map.view`; 
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


