#!/usr/bin/perl

use POSIX;
$L=21.736663658513;
open(LOG, "./md.log") || die "Couldn't find md.log.\n";
while (<LOG>){
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
    if (/Cluster of SiF4/){
	@Retchant="";
	for ($n=0; $n<10; $n++){
	    pop(@others);
	}
	$e=~s/\*|\(|\)//g;
	@e=split(/ /, $e);
	$etchant=$e[5];
	for ($k=0; $k<5; $k++){
	    $a = readline LOG;
	    $a=~s/\*|\(|\)|<|>//g;
	    @a=split(/ /, $a);
	    if ($a[3] == 14){
		$prod=$a[2];
		@Rprod=($a[6], $a[7], $a[8]);
	    }
	    elsif ($a[2]==$etchant){
		@Retchant=($a[6], $a[7], $a[8]);
	    }
	    else{
		push (@others, $a[2]);
	    }
	}
	if (@Retchant==""){
	    print "weird $run\n";
	}
	$nbrs = `cfginfo -tellme $cfg $prod|grep Nbrs`;
	chomp($nbrs);
	$nbrs =~ s/\*|\(|\)|://g;
	@N=split(/ /, $nbrs);
	for ($n=0; $n<=10; $n++){
	    pop(@nbrs);
	}
	if ($N[1]==4){
	    for ($n=0; $n<$N[1]; $n++){
		push (@nbrs, $N[$n+2]);
	    }
	    @seen{@others} = ();
	    foreach $item (@nbrs) {
		if (!exists $seen{$item}){
		    $back = $item;
		}
	    }
	    $c=`cfginfo -tellme $cfg $back|grep \\<`;
	    $c=~s/<|>//g;
	    @c=split(/ /,$c);
	    @Rback=($c[2], $c[3], $c[4]);
	    @Rij=($Rprod[0]-$Rback[0], 
		  $Rprod[1]-$Rback[1], 
		  $Rprod[2]-$Rback[2]);
	    @Rik=($Rprod[0]-$Retchant[0], 
		  $Rprod[1]-$Retchant[1], 
		  $Rprod[2]-$Retchant[2]);
	    if (abs($Rij[0]/$L) > 0.5){
		if ($Rij[0] < 0){ $Rij[0]+=$L; }
		else { $Rij[0]-=$L; }
	    }
	    if (abs($Rij[1]/$L) > 0.5){
		if ($Rij[1] < 0){ $Rij[1]+=$L; }
		else { $Rij[1]-=$L; }
	    }
	    if (abs($Rik[0]/$L) > 0.5){
		if ($Rik[0] < 0){ $Rik[0]+=$L; }
		else { $Rik[0]-=$L; }
	    }
	    if (abs($Rik[1]/$L) > 0.5){
		if ($Rik[1] < 0){ $Rik[1]+=$L; }
		else { $Rik[1]-=$L; }
	    }
	    
	    $Rijmag=sqrt($Rij[0]*$Rij[0]+$Rij[1]*$Rij[1]+$Rij[2]*$Rij[2]);
	    $Rikmag=sqrt($Rik[0]*$Rik[0]+$Rik[1]*$Rik[1]+$Rik[2]*$Rik[2]);
	    if ($Rikmag < 4){
		$dot= $Rij[0]*$Rik[0] + $Rij[1]*$Rik[1] + $Rij[2]*$Rik[2];
		$dot/=($Rijmag*$Rikmag);
		$ang = acos($dot)/3.14159*180;
		print "$run $Rikmag $ang \n";
	    }
	}
    }
}

#open (CFG,$file);
#while (<CFG>){
#    if (/\#|%/){
#	if (/dimensions/){
#	    @L=split;
#	    $Lx=$L[1];
#	    $Ly=$L[2];
#	}
#    }
#    else{
#	s/<|>//g;
#	@atom=split;
#	$ix=$atom[0];
#	$id=$atom[1];
#	$Rx=$atom[2]; $Ry=$atom[3]; $Rz=$atom[4];
#	$Vx=$atom[5]; $Vy=$atom[6]; $Vz=$atom[7];
#	$is_fixed=$atom[8];
#	if ($ix==$etchant){
#	    @Retchant=($Rx,$Ry,$Rz);
#	}
#	if ($ix==$prod){
#	    @Rprod=($Rx,$Ry,$Rz);
#	}
#	if ($ix==$back){
#	    @Rback=($Rx,$Ry,$Rz);
#	}
#    }
#}

