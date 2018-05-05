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
    if (/Cluster/){
	if ($run > $minrun && $run < $maxrun){
	    $ind=index($_,"of ")+3;
	    $cl=substr($_, $ind, index($_,"s")-1-$ind);	
	    if ($cl eq $hsp){
		push @t, $run;
	    }
	}
    }
}
if (!$hist){
    foreach $cl (sort keys %sputprod){
	print "$cl\t$sputprod{$cl}";
	if ($sweepprod{$cl}){
	    print "  ($sweepprod{$cl} swept)";
	}
	print "\n";
    }
}
else{
    for ($i=0; $i<=$#t; $i++){
	print $t[$i], "\n";
    }
}
    
