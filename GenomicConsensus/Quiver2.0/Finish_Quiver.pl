#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt, "quiver:s", "falcon:s", "output:s", "help");


my $help=<<USAGE;
perl $0 --quiver quiver.fasta --falcon ass.fasta
Some of scaffolds do not have reads coverage. These scaffolds will not be polished by Quiver.
This script add these scaffold back to Quiver result.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{output} ||= "falcon_quiver_final.fasta";
my $quiver_seq = getfastaseq_quiver($opt{quiver});
my $falcon_seq = getfastaseq_falcon($opt{falcon});

open OUT, ">$opt{output}" or die "$!";
foreach my $ctg (sort keys %$falcon_seq){
    #print "$c\n";
    if (exists $quiver_seq->{$ctg}){
        my $newseq = formatseq($quiver_seq->{$ctg}, 50);
        #print "quiver\t$ctg\n";
        print OUT ">$ctg\n";
        print OUT "$newseq";
    }else{
        #print "falcon\t$ctg\n";
        my $newseq = formatseq($falcon_seq->{$ctg}, 50);
        print OUT ">$ctg\n";
        print OUT "$newseq";
    }
}
close OUT;

######
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}



sub getfastaseq_falcon
{
$/=">";
my %hash;
my ($file)=@_;
my @chrs;
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
    $hash{$head}=$seq;
    push @chrs, $head;
}
$/="\n";
return \%hash;
}

sub getfastaseq_quiver
{
$/=">";
my %hash;
my ($file)=@_;
my @chrs;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1= $1 if $temp =~/(.*)\|quiver/;
    #print "$temp\t$temp1[0]\n";
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    $hash{$head}=$seq;
    push @chrs, $head;
}
$/="\n";
return \%hash;
}


sub formatseq
{
### format a single line sequence into lines with user specific length
my ($seq, $step)=@_;
my $length=length $seq;
my $run=0;
if ($length % $step == 0){
    $run = int ($length/$step);
}else{
    $run = int ($length/$step) + 1;
}
my $newseq;
for(my $i=0;$i<$run;$i++){
   my $start=$i*$step;
   my $line=substr($seq,$start,$step);
   $newseq.="$line\n";
}
return $newseq;
}


sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}

