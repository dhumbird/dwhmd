#!/usr/bin/perl

$on=0;
open(IN, "./md.log") || die "Couldn't find md.log.\n";
open(OUT, ">md.log.2");
while (<IN>){
    if(/--run /){
	chomp;
	$run=substr($_, index($_ ,"run ")+4);
	$run =~ s/-//g;
	if ($on){ $run+=4000; }
	print OUT "-----------------run $run----------\n";
	if ($run==4000){ $on=1;}
    }
    else{
	print OUT;
    }
}
    
