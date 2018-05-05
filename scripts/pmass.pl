#!/usr/bin/perl

$log=$ARGV[0];
$run=0;

open(LOG, "$log") || die "Couldn't find md.log.\n";
while (<LOG>){
    chomp;
    if(/--run /){
	if ($run!=0) {print "$run $mass $allC\n";}
	$run=substr($_, index($_ ,"run ")+4);
	$run =~ s/-//g;
	$mass=0;
	$allC=0;
    }
    if (/Cluster/){
	$ind=index($_,"of ")+3;
	if (/swept/){
	    $cl=substr($_, $ind, index($_,"swept")-1-$ind);	
	}
	else{
	    $cl=substr($_, $ind, index($_,"sput")-1-$ind);
	}
	$cl=~s/i//;
	$s_ind=index($cl,"S"); 
	$c_ind=index($cl,"C"); 
	$f_ind=index($cl,"F"); 
	$Si=0; $C=0; $F=0;
	if ($cl=~/S/){
	    if ($cl=~/C/){
		$Si=substr($cl,$s_ind+1,$c_ind-1);
	    }
	    elsif ($cl=~/F/){
		$Si=substr($cl,$s_ind+1,$f_ind-1);
	    }
	    else{
		$Si=substr($cl,$s_ind+1,1000);
	    }
	    if ($Si==0) {$Si=1;}
	}
	if ($cl=~/C/){
	    if ($cl=~/F/){
		$C=substr($cl, $c_ind+1, $f_ind-1-$c_ind);
	    }
	    else{
		$C=substr($cl,$c_ind+1,1000);
	    }
	    if ($C==0) {$C=1;}
	}
	if ($cl=~/F/){
	    $F=substr($cl,$f_ind+1,1000);
	    if ($F==0) {$F=1;}
	}
	$allC+=$C;
	$mass+=14*$Si + 6*$C + 9*$F;
    }
}
print "$run $mass $allC\n";
