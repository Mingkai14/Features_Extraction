#!/usr/bin/perl
$libdir=$ARGV[0];
print `./svm-predict ./input_psi.dat $libdir/data/modelpsi_0.005_1.dat ./outputpsi_0.005_1.dat`;
open(fl,"outputpsi_0.005_1.dat");
$i=0;
while($line=<fl>)
{
    $line=~/(\S+)/;
    $score[$i]=$1;
    $i++;
}
close(fl);
$counter=$i;

open(fl,">seq.svr.psi");
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
    if( $i == ($counter-1))
    {
        $tmp=360.0;
    }
    printf fl "%d %8.1f\n",$i+1,$tmp;
}
close(fl);

exit();
