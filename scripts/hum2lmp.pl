#!/usr/bin/perl

$name{14}="Si"; $mass{14}=28.066;
$name{9}="F"; $mass{9}=18.9984;
$name{18}="Ar"; $mass{18}=39.948;
$name{2}="He"; $mass{2}=4.0026;
$name{17}="Cl"; $mass{17}=35.453;
$name{10}="Ne"; $mass{10}=20.17;
$name{6}="C"; $mass{6}=12.011;
$name{1}="H"; $mass{1}=1.0079;

$date=`date`; chomp($date);
$corner = 0;
$buffer = 10;

for($i=0; $i<=$#ARGV; $i++){
	if ($ARGV[$i] =~ /cfg/){
	    $humcfg=$ARGV[$i];
	}
    elsif ($ARGV[$i] =~ /corner/){
		$corner = 1;
    }
}

open(HUM, $humcfg);
while (<HUM>){
	if (!/%/){
		chomp;
		if (/\#/){
		    @data=split;
		    if (/N /){
				$N = $data[1];
		    }
		    if (/Nfixed /){
				$NFix = $data[1];
		    }
		    if (/dimensions/){
				$dx = $data[1];
				$dy = $data[2];
				$dz = $data[3];
		    }
		    if (/time/){
				$t = $data[1];
		    }
		}
	}
}


#$cfacfg="0000.cfg";
$cfacfg=$humcfg;
$cfacfg=~s/cfg/dat/;
open(LMP, ">$cfacfg");

printf LMP "hum2lmp: convert $humcfg to dat at timestep %10.5f\n\n", $t;
print LMP "$N atoms\n";
print LMP "2 atom types\n\n";

printf LMP "%10.5f%10.5f xlo xhi\n",0,$dx;
printf LMP "%10.5f%10.5f ylo yhi\n",0,$dy;
printf LMP "%10.5f%10.5f zlo zhi\n",0,$dz+$buffer;

print LMP "Atoms\n\n";

#print LMP "ITEM: ATOMS id type x y z vx vy vz\n";

$i=1;
open(HUM, $humcfg);
while (<HUM>){
    if (!/%/){
		chomp;
		if (!/\#/){
		    s/<|>//g; 
		    @atom=split;
		    $ix=$atom[0];
		    $id=$atom[1];
		    $Rx=$atom[2]; $Ry=$atom[3]; $Rz=$atom[4];
		    if ($corner==1){
		    	$Rx+=$dx/2+1e-8;
		    	$Ry+=$dy/2+1e-8;
 		    	$Rz+=$dz/2+1e-8;
 		    }
		    $Vx=$atom[5]; $Vy=$atom[6]; $Vz=$atom[7];
		    $is_fixed=$atom[8];
		    print LMP "$i $name{$id} ";
		    printf LMP "%10.5f%10.5f%10.5f",$Rx,$Ry,$Rz;
		    printf LMP "%10.5f%10.5f%10.5f\n",$Vx,$Vy,$Vz;
		    $i++;
	    }
	}
}