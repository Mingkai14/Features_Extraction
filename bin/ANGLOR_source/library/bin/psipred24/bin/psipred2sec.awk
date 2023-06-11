#!/bin/csh -f

#echo "! from $1"
#set base=$1:r
set base=""
cat $1 | awk 'BEGIN{cc="";ir=0}{if(substr($1,1,4)=="Conf") c1=$2;if(substr($1,1,4)=="Pred") {c2=$2;for(i=1;i<61;i++) {d1=substr(c1,i,1);d2=substr(c2,i,1);j=j+1;l=3;if(d2=="H") l=1;if(d2=="E") l=2;if(d1!="") {cc=cc l;ir=ir+1;}}}}END{print ">",ir,"'$base'";for(i=1;i<=ir;i=i+70){print substr(cc,i,70)}}'

