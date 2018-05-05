#!/usr/bin/perl

$launch=1;

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
	$pdb=`/home/dhumbird/bin/pdbdump @ARGV`;
	if ($launch){
	    exec("/usr/local/bin/vmd $pdb")
		or die "Couldn't replace myself with vmd!\n";
	}
    }
    else{
	die("nope\n");
    }
}
else{
    $pdb=`/home/dhumbird/bin/pdbdump @ARGV`;
    if ($launch){
	exec("/usr/local/bin/vmd $pdb")
	    or die "Couldn't replace myself with vmd!\n";
    }
}


    


