#!/usr/bin/perl

if ($#ARGV==0){
    $logfile=$ARGV[0];
}
else{
    $logfile="md.log";
}
open(LOG, $logfile) || die "Couldn't find $logfile.\n"; 

while(<LOG>){
    if (/Cluster/){
	s/\* Atom |sputtered from |\(|\)|<|>|Time: //g;
	push @sput, [ split ];
    }
    
}
#$si_phys=$si_chem=$hal_phys=$hal_chem=0;
#for ($i=0; $i<=$#sput; $i++){
#    $j=$i;
#    while ($i+1<=$#sput && $sput[$i][5]==$sput[$i+1][5]){
#	if ($sput[$i][1]==14){
#	    $si_chem++;
#	}
#	else{
#	    $hal_chem++;
#	}
#	$i++;
#    }
#    if ($j==$i){
#	if ($sput[$i][1]==14){
#	    $si_phys++;
#	}
#	else{
#	    $hal_phys++;
#	}
#    }
#    else{
#	if ($sput[$i][1]==14){
#	    $si_chem++;
#	}
#	else{
#	    $hal_chem++;
#	}
#    }
#}
#print "****** Silicon *****\n";
#print "Chemical:\t$si_chem\n";
#print "Physical:\t$si_phys\n";
#print "****** Halogen *****\n";
#print "Chemical:\t$hal_chem\n";
#print "Physical:\t$hal_phys\n";
#print "--------------------\n";
#print "Total:\t\t",$hal_chem+$hal_phys+$si_chem+$si_phys,"\n";
