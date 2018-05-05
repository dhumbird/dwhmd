#!/usr/bin/perl

for ($i=0; $i<=$#ARGV; $i++){
    if ($ARGV[$i] =~ /cfg/){
	push @files, $ARGV[$i];
    }
} 

for ($i=0; $i<=$#files; $i++){
    print $i,"\t",`cfginfo -amorph $files[$i] |head -1`;
}
