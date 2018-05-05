#!/usr/bin/perl

$name{14}="Si";
$name{9}="F";
$name{18}="Ar";

$humcfg=$ARGV[$#ARGV];
$cfacfg="0000.cfg";
$date=`date`; chomp($date);
open(HUM, $humcfg);
open(CFA, ">$cfacfg");
print CFA "#UnitSystem APVK\n";
print CFA "#Periodicities 1 1 0\n";
$i=0;
while (<HUM>){
    if (!/%/){
	chomp;
	if (/\#/){
	    @data=split;
	    if (/N /){
		print CFA "% Number of atoms in this cfg = $data[1]\n";
	    }
	    if (/Nfixed /){
		print CFA "% Number of Fixed Atoms = $data[1]\n";
	    }
	    if (/dimensions/){
		print CFA "#BoxSize.xyz $data[1] $data[2] $data[3]\n";
	    }
	    if (/time/){
		print CFA "#cfgTime $data[1]\n";
	    }
	}
	else{
	    s/<|>//g;
	    @atom=split;
	    $ix=$atom[0];
	    $id=$atom[1];
	    $Rx=$atom[2]; $Ry=$atom[3]; $Rz=$atom[4];
	    $Vx=$atom[5]; $Vy=$atom[6]; $Vz=$atom[7];
	    $is_fixed=$atom[8];
	    print CFA "$i $name{$id} $Rx $Ry $Rz $Vx $Vy $Vz $is_fixed\n";
	    $i++;
	}
    }
}
