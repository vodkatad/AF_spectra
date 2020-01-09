#!/usr/bin/env perl
use warnings;
use strict;
if (scalar(@ARGV) != 2) {
    die "Usage $0 SNV|indel|both nomulti|multi < headerlessvcf"
}
my $what = $ARGV[0];
my $multi = $ARGV[1];
my $nmulti = 0;
while (<STDIN>) {
    chomp;
    my @line = split("\t", $_);
    die "I expect a single sample vcf" if scalar(@line) != 10;
    my @l = split(':', $line[8]);
    ###FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
    ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
    ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
    die "Wrong vcf FORMAT" if ($l[4] ne 'AD' || $l[3] ne 'RD' || $l[5] ne "FREQ");
    if ($line[4] =~ /,/) {
        if ($multi eq 'multi') {
            die "Sorry still to be implemented"; # probably will need to use a library for this
        }
        $nmulti++;
        next;     
    } else {
        &manage_entry($line[2], $line[3], $line[4], $line[9], $what, $line[0], $line[1]);
    }
}

sub manage_entry {
    my $id = shift;
    my $ref = shift;
    my $alt = shift;
    my $g = shift;
    my $what = shift;
    my $chr = shift;
    my $b = shift;
    $b = $b-1; #switch to zero based
    my $e = $b + length($ref); # end escluded considering length of ref, TODO FIXME for long indels
    #paired varscan:
    #chr1    926278  .       G       T       .       PASS    DP=20;SOMATIC;SS=2;SSC=7;GPV=1;SPV=0.18947;AC=1;AN=2    GT:GQ:DP:RD:AD:FREQ:DP4 0/1:.:9:7:2:22.22%:0,7,1,1
    if ($what eq 'SNV') {
        return if (length($ref) != 1 || length($alt) != 1);
    } elsif ($what eq 'indel') { 
        return if (length($ref) == 1 && length($alt) == 1);
    } # we do not check for both
    my @afs = split(':',$g);
    my $nreadsref = $afs[3];
    my $nreadsalt = $afs[4];
    my $af = $afs[5];
    $af =~ s/%$//;
    $af = $af / 100;
    $id = $id . ':' . $nreadsref . ':' . $nreadsalt . ':' . $af; 
    print $chr . "\t" . $b . "\t" . $e . "\t" .  $id . "\n";
}


print STDERR "multiallelic\t$nmulti";
