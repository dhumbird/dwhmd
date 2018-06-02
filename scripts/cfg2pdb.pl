#!/usr/bin/perl

$launch=1;
$vmd = '/usr/local/bin/vmd';
$pdbdump = '/home/dhumbird/md/bin/pdbdump';

for ($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i]=~/-v/){
	$launch=1;
    }
    elsif ($ARGV[$i]=~/\+v/){
	$launch=0;
    }
}

#print "$#ARGV $ARGV[$#ARGV]\n";

if ($#ARGV==0){
    if (-e $ARGV[$#ARGV]){
	$pdb=`$pdbdump @ARGV`;
	if ($launch){
	    exec("$vmd $pdb")
		or die "Couldn't replace myself with vmd!\n";
	}
    }
    else{
	die("nope\n");
    }
}
else{
    $pdb=`$pdbdump @ARGV`;
    if ($launch){
	exec("$vmd $pdb")
	    or die "Couldn't replace myself with vmd!\n";
    }
}


    


