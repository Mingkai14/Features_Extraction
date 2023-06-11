#!/usr/bin/perl
$path=$ARGV[3];
$ENV{'PATH'}="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/bin/X11:/usr/games:$path";
$ENV{'LD_LIBRARY_PATH'}="/usr/local/lib:/usr/lib:/lib:$path";
print "$path/simple_test $ARGV[0] $ARGV[1] $ARGV[2]\n";
print `$path/simple_test $ARGV[0] $ARGV[1] $ARGV[2]`;
exit();
