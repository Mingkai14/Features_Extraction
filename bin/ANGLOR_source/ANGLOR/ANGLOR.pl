#!/usr/bin/perl
use Math::Trig;

########### setup  the environment and Working DIRectory ###
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/usr/pgi/linux86/bin";
$ENV{'LD_LIBRARY_PATH'}="/usr/local/lib:/usr/lib:/lib";

%ts=(
     'GLY'=>'G',
     'ALA'=>'A',
     'VAL'=>'V',
     'LEU'=>'L',
     'ILE'=>'I',
     'SER'=>'S',
     'THR'=>'T',
     'CYS'=>'C',
     'MET'=>'M',
     'PRO'=>'P',
     'ASP'=>'D',
     'ASN'=>'N',
     'GLU'=>'E',
     'GLN'=>'Q',
     'LYS'=>'K',
     'ARG'=>'R',
     'HIS'=>'H',
     'PHE'=>'F',
     'TYR'=>'Y',
     'TRP'=>'W',

     'ASX'=>'B',
     'GLX'=>'Z',
     'UNK'=>'X',

     'G'=>'GLY',
     'A'=>'ALA',
     'V'=>'VAL',
     'L'=>'LEU',
     'I'=>'ILE',
     'S'=>'SER',
     'T'=>'THR',
     'C'=>'CYS',
     'M'=>'MET',
     'P'=>'PRO',
     'D'=>'ASP',
     'N'=>'ASN',
     'E'=>'GLU',
     'Q'=>'GLN',
     'K'=>'LYS',
     'R'=>'ARG',
     'H'=>'HIS',
     'F'=>'PHE',
     'Y'=>'TYR',
     'W'=>'TRP',

     'a'=>'CYS',
     'b'=>'CYS',
     'c'=>'CYS',
     'd'=>'CYS',
     'e'=>'CYS',
     'f'=>'CYS',
     'g'=>'CYS',
     'h'=>'CYS',
     'i'=>'CYS',
     'j'=>'CYS',
     'k'=>'CYS',
     'l'=>'CYS',
     'm'=>'CYS',
     'n'=>'CYS',
     'o'=>'CYS',
     'p'=>'CYS',
     'q'=>'CYS',
     'r'=>'CYS',
     's'=>'CYS',
     't'=>'CYS',
     'u'=>'CYS',
     'v'=>'CYS',
     'w'=>'CYS',
     'x'=>'CYS',
     'y'=>'CYS',
     'z'=>'CYS',

     'B'=>'ASX',
     'Z'=>'GLX',
     'X'=>'CYS',
    );


$pdb=$ARGV[0]; #pdb name


######Change the followingi settings#######
$datadir="/home/wmk/Python_Pro/Features_Extraction/bin/ANGLOR_source/example/$pdb";
$libdir="/home/wmk/Python_Pro/Features_Extraction/bin/ANGLOR_source/library";
$bindir="$libdir/bin"; #ANGLOR bin directory
$user="wusitao"; # user name
##############End of changes###############################

$blastdir="$bindir/blast/bin"; #Psi-Blast directory
$db="/home/wmk/Python_Pro/Features_Extraction/src/nr_filter/nr.filter";


################# directories #############################
$TAG="ANGLOR_$pdb";
$work_dir="/home/wmk/Python_Pro/Features_Extraction/tmp/";

################ working directory ########################
`/bin/mkdir -p $work_dir`;
chdir "$work_dir";
`/bin/rm -f $work_dir/*`;



################ make fasta sequence file #################
@seqtxts=`cat $datadir/seq.txt`;
$sequence="";
foreach $seqtxt(@seqtxts){
    goto pos6 if($seqtxt=~/\>/);
    $seqtxt=~s/\s//mg;
    $seqtxt=~s/\n//mg;
    $sequence=$sequence.$seqtxt;
  pos6:;
}
$Lch=length $sequence;


open(seq,">protein.seq");
printf seq ">protein\n";
for($i=1;$i<=$Lch;$i++){
    $a=substr($sequence,$i-1,1);
    printf seq "$a";    
    $seqQ3{$i}=$ts{$a};
    
    if($i==int($i/60)*60){
	printf seq "\n";
    }
}
printf seq "\n";
close(seq);

###################################################
#####generate seq.dat :: secondary structure ######
###################################################

#use psipred to predict secondary structure
printf "(1) Runing for secondary structure prediction\n";
print `$libdir/bin/psipred24/runpsipred protein.seq`;


open(psipred,"protein.horiz");
open(yan,">seq.dat");
$j=0;
while($line=<psipred>)
{
    if($line=~/Conf:\s+(\d+)/)
    {
        $conf=$1;
        <psipred>=~/Pred:\s+(\S+)/;
        $pred=$1;
        <psipred>=~/AA:\s+(\S+)/;
        $aa=$1;
        $num=length $aa; 
        for($i=1;$i<=$num;$i++)
        {
            $j++;
            $conf1=substr($conf,$i-1,1);
            $pred1=substr($pred,$i-1,1);
            $aa1=substr($aa,$i-1,1);
            $sec{$j}=1;
            $sec{$j}=2 if($conf1 >=1 && $pred1 eq 'H'); #confidence score >1
            $sec{$j}=4 if($conf1 >=1 && $pred1 eq 'E'); #confidence score >1
            printf yan "%5d %3s %5d %5d\n",$j,$seqQ3{$j},$sec{$j},$conf1;
        }
    }
}
close(yan);
close(psipred);


###############################################################
#make 'exp.dat' :: solvent accessibility
###############################################################
printf "(2) Runing for solvent accessibility prediction\n";
print `cp $libdir/bin/exp.pl ./`;
print `./exp.pl $libdir`;

open(fl,"exp.dat");
<fl>; #skip this line
$in=0;
while($line=<fl>)
{
    $line=~/\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    for($pp=1;$pp<=17;$pp++)
    {
        if($$pp==0)
        {
            last;
        }
    }
    $exp[$in]=($pp-1)/20;#exposed    
    $in++;
}
close(fl);

open(wfl,">seq.exp.sa");
printf wfl "$in\n";
for($i=0;$i<$in;$i++)
{
    printf wfl "%-4d %8.2f\n",$i+1,$exp[$i];
}
close(wfl);


###############################################################
#make '$pdb\_pssm.txt' :: sequence profile
###############################################################
printf "(3) Running Psi-blast .....\n";
`$blastdir/blastpgp  -b 1000 -j 3 -h 0.001 -d $db -i protein.seq -C psitmp.chk -Q $pdb\_pssm.txt > blast.out`;



#####generate seq.svr.psi and seq.svr.phi ######
`cp $libdir/bin/construct_input_torsion.pl ./`;
`cp $libdir/bin/predict_psi.pl ./`;
`cp $libdir/bin/predict_phi_ann.pl ./`;
`cp $libdir/bin/svm-predict ./`;
`cp $libdir/bin/change_format ./`;
`cp $libdir/bin/simple_test.pl ./`;
`cp $libdir/bin/simple_test ./`;

`./construct_input_torsion.pl 10 $pdb`;
`./predict_psi.pl $libdir`;
`./predict_phi_ann.pl $libdir`;


`/bin/cp seq.ann.phi      $datadir/phi.txt`;
`/bin/cp seq.svr.psi      $datadir/psi.txt`;


################# endding procedure ######################
printf "Prediction finished\n";
`sync`;
`sync`;
sleep(1);
`rm -fr $work_dir`;
exit();
