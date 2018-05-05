#!/usr/bin/perl
open(LOG, "md.log");
$save=0;
$i=1;
$a=5.4341659146;
while(<LOG>){
    chomp;
    if (/--run /){
	$run=substr($_, index($_ ,"run ")+4);
	$run =~ s/-//g;
    }
    if (/Adding new/){
	$save=1;
    }
    if (/saved/){
	$cfg=substr($_,index($_,"File: ")+6,1000);
	push (@cfgs, $cfg);
	if ($save==1){
	    $add{$cfg}=$i++;
	    $save=0;
	}
    }
}
close(LOG);

$subt=0;

foreach $file (@cfgs){
    if ($add{$file}!=undef){
	$subt+=$a/2;
    }
    if ($subt!=0){
	if (open(CFG, $file)){
	    $new="tmp".$file;
	    open(NEW, ">$new");
	    while(<CFG>){
		if (/dimensions/){
		    chomp;
		    @dim=split;
		    $Lx=$dim[1];
		    $Ly=$dim[2];
		    $Lz=abs($dim[3]);
		    $Lz-=$subt;
		    print NEW "#dimensions $Lx $Ly $Lz\n";
		}
		else{
		    print NEW $_;
		}
	    }
	    close(CFG);
	    close(NEW);
	    `mv $new $file`;
	}
    }
}
