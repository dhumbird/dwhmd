#!/usr/bin/perl

$ang1=90; $ang2=0; $ang3=0; 
$deg=2;

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
    elsif ($ARGV[$i] =~ /ang/){
	$ang1=$ARGV[++$i];
	$ang2=$ARGV[++$i];
	$ang3=$ARGV[++$i];
    }
    elsif ($ARGV[$i] =~ /deg/){
	$deg=$ARGV[++$i];
    }
    else{
	push (@args, $ARGV[$i]);
    }
}
foreach $file (@files){
    $file=~s/.cfg//;
    `rm -rf $file.rot`;
    `mkdir $file.rot`;
    
    $ang=0;
    while($ang<360){
	$a3=$ang3+$ang;
	`mktiff @args -ang $ang1 $ang2 $a3 $file.cfg`;
	$ANG=$a3;
	while (length($ANG)!=length("360")){
	    $ANG="0".$ANG;
	}
	#`mogrify -geometry 300x300 rot.tif`;
	`mv $file.png $file.rot/$file.$ANG.png`;
	$ang+=$deg;
    }
}





