#!/usr/bin/perl
use POSIX;
if ($#ARGV==0){
    $logfile=$ARGV[0];
}
else{
    $logfile="md.log";
}
open(LOG, $logfile) || die "Couldn't find $logfile.\n"; 
$physical=0;
while(<LOG>){
    chomp;
    if(/--run /){
	$run=substr($_, index($_ ,"run ")+4);
	$run =~ s/-//g;
    }
    if (/Cfg saved/){
	$cfg = substr($_,index($_,"File: ")+6,1000);
    }
    if (/Adding species F./){
	$e=$_;
    }
    if (/Cluster of SiF2 sputtered/){
	for ($n=0; $n<=10; $n++){
	    pop(@others);
	}
	$e=~s/\*|\(|\)//g;
	@e=split(/ /, $e);
	for ($k=0; $k<3; $k++){
	    $a = readline LOG;
	    $a=~s/\*|\(|\)|<|>//g;
	    @a=split(/ /, $a);
	    if ($a[3] == 14){
		$prod=$a[2];
	    }
	    #elsif ($a[2]==$etchant){
	#	@Retchant=($a[6], $a[7], $a[8]);
	    #   }
	    else{
		push (@others, $a[2]);
	    }
	}
	$nbrs = `cfginfo -tellme $cfg $prod|grep Nbrs`;
	chomp($nbrs);
	$nbrs =~ s/\*|\(|\)|://g;
	@N=split(/ /, $nbrs);
	for ($n=0; $n<=10; $n++){
	    pop(@nbrs);
	}
	if ($N[1] > 2){
	    for ($n=0; $n<$N[1]; $n++){
		push (@nbrs, $N[$n+2]);
	    }
	    @seen{@others} = ();
	    $total=0;
	    foreach $item (@nbrs) {
		if ($total >= 0){
		    if (exists $seen{$item}){
			$total++;
		    }
		    else {
			$c=`cfginfo -tellme $cfg $item|grep \\<`;
			$c=~s/<|>//g;
			@c=split(/ /,$c);
			if ($c[1]==9){
			    $total=-10;
			}
		    }
		}
	    }
	    #print "SiF2: $prod @others | @nbrs ";
	    #if ($total==2){ print "physical\n";}
	    #else {print "chemical\n";}
	    if ($total == 2) {$physical++;}
	}
    }
}
print "$physical\n";

#    if (/sput/){
#	s/\* Atom |sputtered from |\(|\)|<|>|Time: //g;
#	push @sput, [ split ];
#    }
#}
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
