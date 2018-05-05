#!/usr/bin/perl

$name{Si}=14;
$name{F}=9;
$name{Ar}=18;
$name{C}=6;

for ($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i]!~/_/){
	push @filelist, $ARGV[$i];
    }
}

#read the origins as the positions of the first file.
#open(CFA, $filelist[0]);
#while(<CFA>){
#    if (!/\#|%/){
#	chomp;
#	@atom=split;
#	$ix=$atom[11];
#	$Ox{$ix}=$atom[2]; $Oy{$ix}=$atom[3]; $Oz{$ix}=$atom[4];
#    }
#}
#close(CFA);

foreach $cfacfg (@filelist){
#    $rhs=$cfacfg=~s/.cfg//;
#    $rhs=$cfacfg%1000;
#    $lhs=($cfacfg-$rhs)/1000;
#    while (length($rhs) < length("999")){
#	$rhs="0".$rhs;
#    }
#    while (length($lhs) < length("9999")){
#	$lhs="0".$lhs;
#    }
#    $humcfg="cfa_".$lhs."-".$rhs.".cfg";
    $humcfg="hum_".$cfacfg;
    $date=`date`; chomp($date);
    open(HUM, ">$humcfg");
    open(CFA, "$cfacfg");

    print HUM "%$date $humcfg\n";
    print HUM "%file format revision 05232002\n";

    while(<CFA>){
	chomp;
	if (/%/){
	    if (/cfg/){
		$N=substr($_,index($_,"=")+1);
	    }
	    if (/Fixed/){
		$Nfixed=substr($_,index($_,"=")+1);
	    }
	}
	if (/#/){
	    if(/BoxSize/){
		@dim=split;
	    }
	}
    }
    print HUM "#N $N\n#Nfixed $Nfixed\n#Nmax $N\n";
    print HUM "#temperature 300.0\n";
    print HUM "#dimensions $dim[1] $dim[2] $dim[3]\n";
    print HUM "#time ", $cfacfg/1000,"\n";
    print HUM "#-----------------\n";
    close(CFA);
    open(CFA, "$cfacfg");
    while(<CFA>){
	if (!/\#|%/){
	    chomp;
	    @atom=split;
	    $ix=$atom[0];
	    $id=$name{$atom[1]};
	    $Rx=$atom[2]; $Ry=$atom[3]; $Rz=$atom[4];
	    $Vx=$atom[5]; $Vy=$atom[6]; $Vz=$atom[7];
	    $is_fixed=$atom[8];
	    #$xpass=$atom[9];
	    #$ypass=$atom[10];
	    print HUM "$ix $id <$Rx $Ry $Rz> <$Vx $Vy $Vz> ";
	    print HUM "$is_fixed\n";
	}
    }
}

