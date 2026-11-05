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
    next if $line[6] ne 'PASS' && $line[6] ne 'recall;alleleBias' && $line[6] ne 'recall;PASS';
    ##recall;PASS
    #PASS
    #recall;alleleBias
    die "I expect a single sample vcf" if scalar(@line) != 10;
    my @l = split(':', $line[8]);
    #die "Wrong vcf FORMAT" if ($line[8] ne 'GT:AD:AF:DP:F1R2:F2R1:SB' && $line[8] ne 'GT:AD:DP:GQ:PL' && $line[8] ne 'GT:AD:AF:F1R2:F2R1:DP:SB:MB');
    #   TODO inefficient split only here...
    die "Wrong vcf FORMAT" if ($l[1] ne 'AD' || $l[2] ne 'AF');
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
    $id = $chr . ":" . $b . ":" . $ref . ":" . $alt;
    #%CHROM:%POS:%REF:%ALT
    $b = $b-1; #switch to zero based
    my $e = $b + length($ref); # end escluded considering length of ref, TODO FIXME for long indels
    ## platypus Giorgio:
    ##chr1    17385   .       G       A       38      recall;MQ;alleleBias    .       GT:AD:AF:DP:F1R2:F2R1:SB        .:38,5:0.116279:43:.:.:.
    if ($what eq 'SNV') {
        return if (length($ref) != 1 || length($alt) != 1);
    } elsif ($what eq 'indel') { 
        return if (length($ref) == 1 && length($alt) == 1);
    } # we do not check for both
    my @afs = split(':',$g);
    my @nreads = split(',',$afs[1]);
    #(snakemake) egrassi@godot:/mnt/trcanmed/snaketree/data/MA_cell_lines_MR/to_MA$ zcat platypus/*.platypus-recall.vcf.gz  |grep -v  "^#" |  cut -f 1,2,3,7,10 | tr ":" "\t" | tr "," "\t"|bawk '$7==0'
    # probably not required, but...
    if ($nreads[1] != 0) {
        $id = $id . ':' . $nreads[0] . ':' . $nreads[1] . ':' . $afs[2]; 
        print $chr . "\t" . $b . "\t" . $e . "\t" .  $id . "\n";
    }
}


print STDERR "multiallelic\t$nmulti";
