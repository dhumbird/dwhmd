#!/usr/bin/perl

$minrun=0;
$maxrun=1e18;
$hist=0;
for ($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i]=~/max/){
	$maxrun=$ARGV[++$i];
    }
    if ($ARGV[$i]=~/min/){
	$minrun=$ARGV[++$i];
    }
    if ($ARGV[$i]=~/hist/){
	$hist=1;
	$hsp=$ARGV[++$i];
    }
}
open(LOG, "./md.log") || die "Couldn't find md.log.\n";
while (<LOG>){
    chomp;
    if(/--run /){
	$run=substr($_, index($_ ,"run ")+4);
	$run =~ s/-//g;
    }
    if (/Adding species/){
	if (/F|Cl/){
	    $thermal=1;
	}
	if (/Ar/){
	    $thermal=0;
	}
    }
    if (/Cluster/){
	if ($run > $minrun && $run < $maxrun){
	    $ind=index($_,"of ")+3;
	    if (/swept/){
		$cl=substr($_, $ind, index($_,"swept")-1-$ind);	
		$sweepprod{$cl}++;
	    }
	    else{
		$cl=substr($_, $ind, index($_,"sput")-1-$ind);
	    }
	    $sputprod{$cl}++;
	    if ($cl=~/$hsp/){
		push @t, $run;
	    }
	    if ($thermal){
		$thermprod{$cl}++;
	    }
	}
    }
}
if (!$hist){
    foreach $cl (sort keys %sputprod){
	print "$cl\t$sputprod{$cl}";
	if ($sweepprod{$cl}){
	    print " $sweepprod{$cl} swept";
	}
	if ($thermprod{$cl}){
	    print " $thermprod{$cl} thermal";
	}
	print "\n";
    }
}
else{
    for ($i=0; $i<=$#t; $i++){
	print $t[$i], "\n";
    }
}
    
