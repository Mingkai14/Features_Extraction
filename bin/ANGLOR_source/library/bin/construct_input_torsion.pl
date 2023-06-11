#!/usr/bin/perl 
%SEQ=(
        "ALA"=>"A",
        "ARG"=>"R",
        "ASN"=>"N",
        "ASP"=>"D",
        "CYS"=>"C",
        "GLN"=>"Q",
        "GLU"=>"E",
        "GLY"=>"G",
        "HIS"=>"H",
        "ILE"=>"I",
        "LEU"=>"L",
        "LYS"=>"K",
        "MET"=>"M",
        "PHE"=>"F",
        "PRO"=>"P",
        "SER"=>"S",
        "THR"=>"T",
        "TRP"=>"W",
        "TYR"=>"Y",
        "VAL"=>"V",
        );
$window=$ARGV[0]; #  10
$pdb{'0'}=$ARGV[1]; # pdb name
$counter=1;

open(fl3,">input_phi.dat");
open(fl4,">input_psi.dat");
for($i=0;$i<$counter;$i++)
{
    printf "$i :: $pdb{$i}\n";
##########################################
#read from exp.dat, from 0,1,2,3,...
###########################################
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
        $exp{$i,$in}=($pp-1)/20;#exposed
        
	if($exp{$i,$in}<0.25)#buried
        {
	    $exposed{$i,$in}=0;
	    $buried{$i,$in}=1;
	}
	else #exposed
	{
	    $exposed{$i,$in}=1;
	    $buried{$i,$in}=0;
	}       
        $in++;

    }
    close(fl);
#################

################################
#read from *_pssm.txt,  from 0,1,2,3,...
################################
    print "$pdb{$i}\_pssm.txt\n";
    open(fl,"$pdb{$i}\_pssm.txt");
    <fl>;#skip this line
    <fl>;#skip this line
    <fl>;#skip this line
    $j=0;

    while(1)
    {
        $line=<fl>;
        $temp_str=substr($line,11,1);        
        $line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
        last if($temp_str eq "");
        $pos{$i,$j}=$1;
        $seq{$i,$j}=$2;       
        for($k=0;$k<20;$k++)
        {
            $kk=$k+3;
            $pssm{$i,$j,$k}="$$kk";
            $pssm{$i,$j,$k}/=10.0;
            $pssm{$i,$j,$k}+=1.0;
            $pssm{$i,$j,$k}/=2.0;            

        }      
        $j++;
    }
    $len_protein=$j;
    close(fl);

 

##############################################################
#read from seq.dat (secondary structure) from 0,1,2,3,...
##############################################################
    print "seq.dat\n";
    open(fl,"seq.dat");
    $j=0;
    while($line=<fl>)
    {
	$line=~/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	$seq_list{$i,$j}=$SEQ{$2};
	print "residue name is different!\n" if($seq_list{$i,$j} ne $seq{$i,$j});
	for($pp=1;$pp<=3;$pp++)
	{
	    $ss{$i,$j,$pp}=0;#coil
	}
	if($3==1)
	{
	    $ss{$i,$j,1}=1;#coil
            $ss2{$i,$j}=0;
	}
	elsif($3==2)
	{
	    $ss{$i,$j,2}=1;#helix
            $ss2{$i,$j}=1;
	}
	elsif($3==4)
	{
	    $ss{$i,$j,3}=1;#beta
            $ss2{$i,$j}=-1;
	}
	$j++;
    }
    close(fl);
    print "protein length is different!\n" if($len_protein!=$j);




####################################
#generate input for training
####################################

    for($j=0;$j<$len_protein;$j++)
    {
	printf fl3 "%-8.3f ", 0;
	printf fl4 "%-8.3f ", 0;
	$cc=0;
	for($k=-$window;$k<=$window;$k++)
	{
	    $kk=$k+$j;            
	    if($kk<0 || $kk>=$len_protein)#blank residue
	    {
		for($p=1;$p<=3;$p++)
		{
		    $cc++;		    
		}
		for($p=0;$p<20;$p++)
		{
		    $cc++;		   
		}
		$cc++;	
                $cc++;                
	    }
	    else
	    {
		for($p=1;$p<=3;$p++)
                {
                    $cc++;
                    printf fl3 "$cc:%1d ",$ss{$i,$kk,$p} if($ss{$i,$kk,$p}!=0);
		    printf fl4 "$cc:%1d ",$ss{$i,$kk,$p} if($ss{$i,$kk,$p}!=0);
                }
                for($p=0;$p<20;$p++)
                {
                    $cc++;
                    printf fl3 "$cc:%-5.3f ", $pssm{$i,$kk,$p} if($pssm{$i,$kk,$p}!=0);
		    printf fl4 "$cc:%-5.3f ", $pssm{$i,$kk,$p} if($pssm{$i,$kk,$p}!=0);                                        
                }                
		$cc++;
		printf fl3 "$cc:%1d ",$exposed{$i,$kk} if($exposed{$i,$kk}!=0);
		printf fl4 "$cc:%1d ",$exposed{$i,$kk} if($exposed{$i,$kk}!=0);
                $cc++;
                printf fl3 "$cc:%1d ",$buried{$i,$kk} if($buried{$i,$kk}!=0);
		printf fl4 "$cc:%1d ",$buried{$i,$kk} if($buried{$i,$kk}!=0);
	    }
	}
	printf fl3 "\n";
	printf fl4 "\n";
    }
    

}#end i
close(fl3);
close(fl4);
exit();
