#!/usr/bin/env perl
use warnings;
use strict;
if (scalar(@ARGV) != 0) {
    die "Usage $0 < pseudompileup"
}
my $tot = 0;
my $totmut = 0;
my $chr = '';
my $b = 0;
my $e = 0;
my $laste = 0;
my $lastchr = '';
my $already = 0;
while (<STDIN>) {
    chomp;
    my @line = split("\t", $_);
    my $bases = $line[5];
    # we must manage end of chr if we ended it in a non mutated interval
    # but the management needs to be done depending on how we begin the new chr. ?
    if ($lastchr ne '' && $line[0] ne $lastchr && $b != 0) {
        print $lastchr . "\t" . $b . "\t" . $laste . "\n";
        $already = 1;
        $b = 0;
    }
    if ($bases =~ /[ACGTNacgtn><*]/) {
        $totmut++;
        #end of a non mutated interval (or going on with muts)
        if ($b != 0) {
            if ($already == 0) {
                $e = $line[1];
                $chr = $line[0]; 
                print $chr . "\t" . $b . "\t" . $e . "\n";
            }
            $chr = '';
            $b = 0;
            $e = 0;
        }
    } else {
        if ($b == 0) { # we open a non mutated interval
            $chr = $line[0];
            $b = $line[1];
        } 
        # if $b !=0 we are extending our non mutated interval
    }
    $already = 0;
    $tot++;
    $laste = $line[2];
    $lastchr = $line[0];
}
# last interval, we assume it will be on the last seen chr.
if ($b != 0) {
    print $lastchr . "\t" . $b . "\t" . $laste . "\n";
}
my $totwt = $tot - $totmut;
print STDERR "$tot\t$totmut\t$totwt\n";