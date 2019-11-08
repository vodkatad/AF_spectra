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
    die "Wrong vcf FORMAT" if ($line[8] ne 'GT:AD:AF:DP:F1R2:F2R1:SB' && $line[8] ne 'GT:AD:DP:GQ:PL');
    if ($line[4] =~ /,/) {
        if ($multi eq 'multi') {
            die "Sorry still to be implemented"; # probably will need to use a library for this
        }
        $nmulti++;
        next;     
    } else {
        &manage_entry($line[2], $line[3], $line[4], $line[9], $what);
    }
}

sub manage_entry {
    my $id = shift;
    my $ref = shift;
    my $alt = shift;
    my $g = shift;
    my $what = shift;
    #mutect:
    #GT:AD:AF:DP:F1R2:F2R1   0/1:14,4:0.235:18:7,2:7,2 
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="Count of reads in F1R2 pair orientation supporting each allele">
    ##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="Count of reads in F2R1 pair orientation supporting each allele">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #haplotypecaller:
    #GT:AD:DP:GQ:PL  0/1:365,53:486:99:132,0,7982
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    # both good
    if ($what eq 'SNV') {
        return if (length($ref) != 1 || length($alt) != 1);
    } elsif ($what eq 'indel') { 
        return if (length($ref) == 1 && length($alt) == 1);
    } # we do not check for both
    my @afs = split(":",$g);
    my @nreads = split(',',$afs[1]);
    print $id . "\t" . $afs[0] . "\t" . $nreads[0] . "\t" . $nreads[1] . "\t" . $afs[2] . "\n";
}


print STDERR "multiallelic\t$nmulti";