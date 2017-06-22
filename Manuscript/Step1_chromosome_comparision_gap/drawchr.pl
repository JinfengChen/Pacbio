#!/usr/bin/perl

use strict;
use FindBin qw ($Bin);
use Getopt::Long;
use SVG;

#################################Usage ##############################################
my %opt;
GetOptions(\%opt,"align:s", "agp:s","help:s");

my $help=<<USAGE;
perl $0 --align --agp

USAGE

if(defined $opt{help} or keys %opt < 1){
        die  $help ;
}

my $qrychr=readtable("CLM.chr.chrlen");
my $qrygap=readagp("CLM.chr.agp");
#my $qrycent=readtable("../input/gramenev1.4.cent");
my $refchr=readtable("FCM.chrlen");
my $refgap=readagp("FCM.chr.agp");
#my $refcent=readtable("../input/tigr6.cent");
my $align=readalign($opt{align});

my $svg=SVG->new(width=>1600, height=>1900);
my $startw=130; 
my $endw  =1500;
my $starth=130;
my $endh  =1700;

my $scale=getscale($refchr,"1400");


#### chromosome line
my %refpos;
my %qrypos;
my @chr=sort keys %$refchr;
for (my $i=0;$i<@chr;$i++){
    ###################reference chromosome
    my $refy=$i*200+$starth;
    my $refx=$startw;
    my $refh=20;
    my $refw=$refchr->{$chr[$i]}/$scale;
    $refpos{$chr[$i]}=$refy;
    #print "$i\t$chr[$i]\t$refx\n";
    my $refline=$svg->rectangle(
              x=>$refx,y=>$refy,
              width=>$refw,height=>$refh,
              rx=>15,ry=>15,
              style=>{
                   fill=>'none', stroke=>'red'
                   #stroke=>'#7A8B8B'
              }
              
    );
    #gap
    my $gap_color="black";
    foreach my $gap (@{$refgap->{$chr[$i]}}){ 
        my $refgapx=$refx + $gap->[0]/$scale;
        my $refgapy=$refy;
        my $refgapw=($gap->[1]-$gap->[0])/$scale;     
        my $refgaph=20;
        my $refgapbox=$svg->rectangle(
                x=>$refgapx,y=>$refgapy,
                width=>$refgapw,height=>$refgaph,
                style=>{
                     fill=>$gap_color, stroke=>$gap_color, 'stroke-width'=>0.1
                }
        );
    }

    #my $refcx=$refx+$refw/2;
    #my $refcy=$refcent->{$chr[$i]}/$scale+$refy;
    #my $refrx=$refw/2;
    #my $refry=$refrx;
    #$svg=centromere($svg,$refcx,$refcy,$refrx,$refry,"black");
    #################query chromosome
    my $qryy=$refy+110;
    my $qryx=$refx;
    my $qryh=$refh;
    my $qryw=$qrychr->{$chr[$i]}/$scale;
    $qrypos{$chr[$i]}=$qryy;
    my $qryline=$svg->rectangle(
              x=>$qryx,y=>$qryy,
              width=>$qryw,height=>$qryh,
              rx=>15,ry=>15,
              style=>{
                   fill=>'none', stroke=>'blue'
              }
    );
    #gap
    foreach my $gap (@{$qrygap->{$chr[$i]}}){ 
        my $qrygapx=$qryx + $gap->[0]/$scale;
        my $qrygapy=$qryy;
        my $qrygapw=10*($gap->[1]-$gap->[0])/$scale;
        my $qrygaph=20;
        my $qrygapbox=$svg->rectangle(
                x=>$qrygapx,y=>$qrygapy,
                width=>$qrygapw,height=>$qrygaph,
                style=>{
                     fill=>$gap_color, stroke=>$gap_color, 'stroke-width'=>0
                       }
             );
    }

    #my $qrycx=$qryx+$qryw/2;
    #my $qrycy=$qrycent->{$chr[$i]}/$scale+$qryy;
    #my $qryrx=$qryw/2;
    #my $qryry=$qryrx;
    #$svg=centromere($svg,$qrycx,$qrycy,$qryrx,$qryry,"black");

    #name 
    my $refnote =$svg->text(
            x=>30, y=>$refy+80,
                style=>{
                     'font-size'=>'32','text-anchor'=>'start','stroke-width'=>0.1
                }
       )->cdata("$chr[$i]"); 

 
    #draw align
    my $color = 'gray';
    foreach my $hit (@{$align->{$chr[$i]}}){
        my $match_up_x1 =$startw + $hit->[0]/$scale;
        my $match_up_x2 =$startw + $hit->[1]/$scale;
        my $match_down_x1 = $startw + $hit->[2]/$scale;
        my $match_down_x2 = $startw + $hit->[3]/$scale;
        my $y1 = $refy + $refh + 5;
        my $y2 = $qryy - 5;
        my $xv = [$match_up_x1,$match_up_x2,$match_down_x2,$match_down_x1];
        my $yv = [$y1, $y1, $y2, $y2];
        my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
        my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color,
                        'fill-opacity' => 0.7
                     }
              );
    }

}


=pod
my $refstarth=$starth+200;
my $refwidth =$endw-$startw;
my $refline  =$svg->rectangle(
              x=>$startw,y=>$refstarth,
              width=>$refwidth,height=>10,
              rx=>5.2,ry=>2.4,
              style=>{
                  stroke=>'black'
              }
);
my $refnote  =$svg->text(
              x=>60, y=>$refstarth+8,
              style=>{
                   fontsize=>'2','text-anchor'=>'start','stroke-width'=>0.1
              }
)->cdata("test");
=cut
writesvg("test.svg",$svg);

#####

sub getscale
{
my ($chr,$height)=@_;
my @len=sort {$b <=> $a} values %$chr;
print "$len[0]\n";
my $rate=$len[0]/$height;
return $rate;
}


sub centromere
{
my ($svg,$cx,$cy,$rx,$ry,$fillcolor)=@_;
my $tag = $svg->ellipse(
        cx=>$cx, cy=>$cy,
        rx=>$rx, ry=>$ry,
        style=>{
            'fill'=>$fillcolor,
        }
    );
return $svg;
}


sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

#  96076   112046  |    55348    71255  |    15971    15908  |    97.65  | chr1	scaffold_1
#    115346   126384  |    73734    84741  |    11039    11008  |    98.17  | chr1	scaffold_1
sub readalign
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my $line = $_;
    $line =~ s/^\s+//;
    $line =~ s/\s+/\t/g;
    my @unit=split("\t", $line);
    #print "$line\n$unit[0]\t$unit[1]\t$unit[3]\t$unit[4]\n";
    #print "$unit[11]\t$unit[0]\t$unit[1]\t$unit[3]\t$unit[4]\n";
    push (@{$hash{$unit[11]}}, [$unit[0], $unit[1], $unit[3], $unit[4]]);
}
close IN;
return \%hash;
}

## Generated from Velvet assembly file /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.chr.fa using script /rhome/cjinfeng/BigData/software/bin/
#chr1	1	49224	0	W	contig_1	1	49224	+
#chr1	49225	49324	1	N	100	fragment	yes	
sub readagp
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my $line = $_;
    my @unit=split("\t", $line);
    if ($unit[4] eq "N"){
       push (@{$hash{$unit[0]}}, [$unit[1], $unit[2]]);
    }
}
return \%hash;
}

################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
system("/rhome/cjinfeng/BigData/software/draw/svg2xxx_release/svg2xxx $file -t pdf");
}


