#!/usr/bin/perl

$XNT1=$ARGV[0];
$XNT2=$ARGV[1];
$CONJUG=$ARGV[2];

$ic=0;
for ($i=1; $i<=4; $i++){
    for ($j=1; $j<=4; $j++){
	for ($k=1; $k<=4; $k++){
	    $ic++;
	    $in3[$ic][1]=$i-1;
	    $in3[$ic][2]=$j-1;
	    $in3[$ic][3]=$k-1;
	}
    }
}
#open(CC, "< /home/humbird/MD/Spline/inter3dtors.d") || die;  $var="tcc";
open(CC, "< /home/humbird/MD/Spline/inter3d_iv_new.d") || die; $var="fcc";
#open(CC, "< /home/humbird/MD/Spline/inter3d_ch.d") || die; $var="fch";
#open(CC, "< /home/humbird/MD/Spline/inter3d_h.d") || die; $var="fhh";

readline *CC;

for ($p=1; $p<=144; $p++){
    $line=readline *CC;
    chomp($line);
    @LMN=split(/\s+/, $line);
    $L=$LMN[1]; 
    $M=$LMN[2]; 
    $N=$LMN[3];
    for ($I=1; $I<=63; $I+=3){
	$line=readline *CC;
	chomp($line);
	@LMN=split(/\s+/, $line);
	$CLMN[$L][$M][$N][$I]=$LMN[1];
	$CLMN[$L][$M][$N][$I+1]=$LMN[2];
	$CLMN[$L][$M][$N][$I+2]=$LMN[3];
    }
    $line=readline *CC;
    chomp($line);
    $line=~s/^\s+//;
    $line=~s/\s+$//;
    $CLMN[$L][$M][$N][64]=$line;
}

for ($L=1; $L<=10; $L++){
    for ($M=1; $M<=10; $M++){
	for ($I=1; $I<=64; $I++){
	    $CLMN[$L][$M][10][$I]=$CLMN[$L][$M][9][$I];
	}
    }
}

#$XNT1=2; $XNT2=3; $CONJUG=5; # my units

$XNT1+=1; $XNT2+=1; # brenner's
$L=int($XNT1); if ($L>=4){$L=4;$XNT1=4;}
$M=int($XNT2); if ($M>=4){$M=4;$XNT2=4;}
$N=int($CONJUG); if ($N>=9){$N=9;$CONJUG=9;}
$tcc=0;
$tcc1=0;
$tcc2=0;
$tcc3=0;
for ($J=1; $J<=64; $J++){
    $x=$CLMN[$L][$M][$N][$J]*($XNT1**$in3[$J][1])
	*($XNT2**$in3[$J][2])
	*($CONJUG**$in3[$J][3]);
    $tcc+=$x;
    $tcc1+=$x*$in3[$J][1]/$XNT1;
    $tcc2+=$x*$in3[$J][2]/$XNT2;
    $tcc3+=$x*$in3[$J][3]/$CONJUG;
}
if (abs($tcc)<1e-7){$tcc=0;}
if (abs($tcc1)<1e-7){$tcc1=0;}
if (abs($tcc2)<1e-7){$tcc2=0;}
if (abs($tcc3)<1e-7){$tcc3=0;}

print "$tcc $tcc1 $tcc2 $tcc3\n";

