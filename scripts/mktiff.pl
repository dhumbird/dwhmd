#!/usr/bin/perl

$cfg2r3d="/home/humbird/md/bin/cfg2r3d";
#$cfg2r3d="/home/humbird/md/unstable/g/cfg2r3d";
$rotate=0;
$labels=0;
$every=1;
$convert="";
$amorph=0;
@ang=(90,0,0);
$tmpcfg="sctemp_00000000.cfg";
for($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i] =~ /\b\d{8}\b/){
	$file=`ls *_$ARGV[$i].cfg* 2>/dev/null`;
	chomp($file);
	if ($file =~ /cfg/){
	    push (@files, $file);
	}
	else{
	    push (@args, $ARGV[$i]);
	}
    }
    elsif ($ARGV[$i] =~ /cfg/){
	push (@files, $ARGV[$i]);
    }
    elsif ($ARGV[$i] =~ /rotate/){
	$rotate=$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /every/){
	$every=$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /convert/){
	$convert=$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /amorph/){
	$amorph=1;
    }
    elsif ($ARGV[$i] =~ /label/){
	$labels=1;
	$labelfile=$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-ang/){
	$ang[0] = $ARGV[++$i];
	$ang[1] = $ARGV[++$i];
	$ang[2] = $ARGV[++$i];
    }
    else{
	push (@args, $ARGV[$i]);
    }
}

for($i=0; $i<=$#files; $i++){
    $ang="";
    $file=$files[$i];
    if ($i % $every == 0){
	if ($amorph){
	    push (@args, "-atri");
	    $atri=`cfginfo -amorph $file -notrans -cutoff 0.35`;
	    chomp $atri;
	    push (@args, $atri);
	}
	$file=substr($file,0,-4);
	if ($rotate!=0){
	    $ang[2]=$rotate*($i/$#files);
	}
	`$cfg2r3d @args -ang @ang $file.cfg >temp.r3d`;
	if ($labels){
	    `cat temp.r3d $labelfile |label3d -png $file.png 2>/dev/null`;
	}
	else{
	    `cat temp.r3d | render -png $file.png 2>/dev/null`;
	}
	`rm -f temp.r3d`;
	if ($convert ne ""){
	    `convert $file.png $file.$convert`;
	    `rm $file.png`;
	}
    }
}
