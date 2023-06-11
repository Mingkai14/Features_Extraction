#!/usr/bin/perl
$libdir=$ARGV[0];

print `./change_format input_phi.dat input_phi_ann.dat 525 1`;
print `./simple_test.pl input_phi_ann.dat $libdir/data/model_ann_angle_phi.dat_50_1000 output_phi_ann.dat $libdir/bin\n`;


open(fl,"output_phi_ann.dat");
$i=0;
while($line=<fl>)
{
    $line=~/\S+\s+(\S+)/;
    $score[$i]=$1;
    $i++;
}
close(fl);
$counter=$i;

open(fl,">seq.ann.phi");
for($i=0;$i<$counter;$i++)
{   
    $tmp=$score[$i]*360;
    if($tmp>180)
    {
	$tmp-=360;
    }
    elsif($tmp<-180)
    {
	$tmp+=360;
    }
    if( $i == 0)
    {
        $tmp=360.0;
    }
    printf fl "%d %8.1f\n",$i+1,$tmp;
}
close(fl);
exit();



