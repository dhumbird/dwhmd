#!/usr/bin/perl

$name{14}="Si"; 
$name{9}="F"; 
$name{18}="Ar"; 
$name{2}="He"; 
$name{17}="Cl";
$name{10}="Ne"; 
$name{6}="C";
$name{1}="H"; 

$i=0;
while($i<=$#ARGV){
    $infile=@ARGV[$i];
    $outfile=$infile.2;
    open(IN, $infile) || die "couldn't open file\n";
    open(OUT, ">$outfile") || die "couldn't open file\n";
    $date=`date`; chomp($date);
    print OUT "%$date $infile\n";
    print OUT "%file format revision 05232002\n";
    while(<IN>){
	if (/file/){
	    s/\D//g;
	    $rev=$_;
	}
    }
    close(IN);
    open(IN, $infile);
    if ($rev==undef){
	while(<IN>){
	    if (!/%/){
		s/Norig/Nmax/;
		print OUT $_; 
	    }
	}
    }
    if ($rev==0){
	while(<IN>){
	    if (!/%/){
		s/Norig/Nmax/;
		print OUT $_; 
	    }
	}
    }
    if ($rev="04152001"){
	while(<IN>){
	    if (!/%/){
		if (/\#/){
		    if (!/material/){
			print OUT $_;
		    }
		}
		else{
		    s/<|>//g;
		    @atom=split;
		    $ix=$atom[0];
		    $id=$atom[1];
		    $Rx=$atom[2]; $Ry=$atom[3]; $Rz=$atom[4];
		    $Vx=$atom[5]; $Vy=$atom[6]; $Vz=$atom[7];
		    $is_fixed=$atom[11];
		    print OUT "$ix $id <$Rx $Ry $Rz> <$Vx $Vy $Vz> $is_fixed\n";
		}
	    }
	}
    }
    `mv -f $outfile $infile`;
    close(IN);
    close(OUT);
    $i++;
}
