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
    #die "Wrong vcf FORMAT" if ($line[8] ne 'GT:AD:AF:DP:F1R2:F2R1:SB' && $line[8] ne 'GT:AD:DP:GQ:PL' && $line[8] ne 'GT:AD:AF:F1R2:F2R1:DP:SB:MB');
    #   TODO inefficient split only here...
    die "Wrong vcf FORMAT" if ($l[4] ne 'NR' || $l[5] ne 'NV');
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
    if ($what eq 'SNV') {
        return if (length($ref) != 1 || length($alt) != 1);
    } elsif ($what eq 'indel') { 
        return if (length($ref) == 1 && length($alt) == 1);
    } # we do not check for both
    my @afs = split(':',$g);
    my $af = 0;
    if ($afs[4] != 0) {
        $af = $afs[5]/$afs[4];
    }
    if ($afs[5] != 0) {
        $id = $id .'@'.$afs[0] . ':' . $afs[4] . ':' . $afs[5] . ':' . $af; 
        print $chr . "\t" . $b . "\t" . $e . "\t" .  $id . "\n";
    }
}


print STDERR "multiallelic\t$nmulti";
