#!/usr/bin/perl

$tmpcfg="sctemp_00000000.cfg";
$stereo=0;
$cfg2r3d="/home/humbird/md/bin/cfg2r3d";
$labels=0;

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
    elsif ($ARGV[$i] =~ /label/){
	$labels=1;
	$labelfile=$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /-ste/) {$stereo=1;}
    else{
	push (@args, $ARGV[$i]);
    }
}

if (!$stereo){
    foreach $file (@files){
	$file=substr($file,0,-4);
	`$cfg2r3d @args $file.cfg >temp.r3d`;
	if ($labels){
	    `cat temp.r3d $labelfile |label3d -png $file.png 2>/dev/null`;
	}
	else{
	    `cat temp.r3d | render -png $file.png 2>/dev/null`;
	}
	`rm -f temp.r3d`;
	exec("display $file.png &");
    }
}

else{
    foreach $file (@files){
	`cfg2r3d @args $file >showcfg.r3d`;
	`stereo3d showcfg.r3d 2>/dev/null`;
	`rm -f showcfg.r3d`;
	exec("display stereo.tiff &");
    }
}
