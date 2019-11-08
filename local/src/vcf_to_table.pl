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
    my @line = $_;
    die "I expect a single sample vcf" if scalar(@line) != 10;
    die "Wrong vcf FORMAT" if $line[0] ne "GT:AD:AF:DP:F1R2:F2R1";
    if ($line[4] =~ /,/) {
        if ($multi eq 'multi') {
            die "Sorry still to be implemented"; # probably will need to use a library for this
        }
        $nmulti++;
        next;     
    } else {
        &manage_nomulti($line[2], $line[3], $line[4], $line[9], $what);
    }
}

sub manage_entry($id, $ref, $alt, $g, $what) {
    #GT:AD:AF:DP:F1R2:F2R1   0/1:14,4:0.235:18:7,2:7,2 
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele">
    ##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    if ($what eq 'SNV') {
        return if (length($ref) != 1 || length($alt) != 1);
    } elsif ($what eq 'SNV')
        return if (length($ref) == 1 && length($alt) == 1);
    } # we do not check for both
    @afs = split(":",$g);
    my @nreads = split(',',$afs[1]);
    print $id . "\t" . $afs[0] . "\t" . $nreads[0] . "\t" . $nreads[1] . "\t" . $afs[1] . "\n";
}


print STDERR "multiallelic\t$nmulti";