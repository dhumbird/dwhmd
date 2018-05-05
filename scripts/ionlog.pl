#!/usr/bin/perl

$cfgfile="temp_00000000.cfg";
$res=1.4;

for ($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i] =~ /cfg/) {
	$cfgfile=$ARGV[$i];
    }
    elsif ($ARGV[$i] =~ /res/){
	$res=$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /log/){
	push (@logs, $ARGV[$i]);
    }
}
 
open(CFG, $cfgfile);
while (<CFG>){
    if (/dimensions/){
	@dim=split;
	$Lx=$dim[1];
	$Ly=$dim[2];
	$Lz=abs($dim[3]);
	$Lx2=$Lx/2;
	$Ly2=$Ly/2;
	$Lz2=$Lz/2;
    }
#    if (!/#|%/){
#	s/<|>//g;
#	@atom=split;
#	$Rx=$atom[2];
#	$Rz=$atom[4];
#	$i=int(($Rx+$Lx2)/$res);
#	if ( ($Lz2-$Rz < $profile[$i]) || $profile[$i]==undef){
#	    $profile[$i]=$Lz2-$Rz;
#	}
#    }
}

$Nx=int($Lx/$res)+1;
$Ny=int($Ly/$res)+1;
for ($i=0; $i<$Nx; $i++){
    for ($j=0; $j<$Ny; $j++){
	$ions[$i][$j]=0;
    }
}

foreach $logfile (@logs){
    open(LOG, $logfile);
    while(<LOG>){
	chomp;
	if (/--run /){
	    $run=substr($_, index($_ ,"run ")+4);
	    $run =~ s/-//g;
	}
	if (/Species added/){
	    $_=substr($_,index($_,"<")+1,index($_,">")-index($_,"<")-1);
	    ($ionx, $iony, $ionz)=split;
	    $ions[int(($ionx+$Lx2)/$res)][int(($iony+$Ly2)/$res)]++;
	}
#	if (/sputtered/){
#	    $_=substr($_,index($_,"<")+1,index($_,">")-index($_,"<")-1);
#	    ($sputx, $sputy, $sputz)=split;
#	    $sputs[int(($sputx+$Lx2)/$res)]++;
#	}
    }
    close(LOG);
}

for ($i=0; $i<$Nx-1; $i++){
    for ($j=0; $j<$Ny-1; $j++){
	if ($ions[$i][$j]==undef){
	    printf "%10.4f\t%10.4f\t", $res*$i-$Lx2, $res*$j-$Ly2;
	    printf "%10.4f\n",0;
	}
	else{
	    printf "%10.4f\t%10.4f\t", $res*$i-$Lx2, $res*$j-$Ly2;
	    printf "%10.4f\n", $ions[$i][$j];
	}
    }
}

