#!/usr/bin/perl

$cfg1=$ARGV[0];
$cfg2=$ARGV[1];
$_=`fcmd -cgmin -nosave $cfg1 2>/dev/null |tail -1`;
@a=split;
$e1=@a[1];
$_=`fcmd -cgmin -nosave $cfg2 2>/dev/null |tail -1`;
@a=split;
$e2=@a[1];
print $e2-$e1 ,"\n";
