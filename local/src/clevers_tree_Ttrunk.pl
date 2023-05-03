#!/usr/bin/env perl

# the logic here was simpler and following clevers trees, each mut has an annoted MRCA and we keep only those whose MRCA is...
#Nodes of the phylogenetic trees (indicated in the column “mrca”) correspond to the labels in the trees in extended data figure 4. For the truncal mutations of P2, substitutions have been assigned prior to the WGD (“17_prewgd”), posterior to the WGD(“17_postwgd”) or unspecified (“17_not_timed”)

#https://www.nature.com/articles/s41586-018-0024-3#Sec44

#https://www.nature.com/articles/s41586-018-0024-3/figures/7
use File::Copy;
use warnings;
use strict;

if (scalar(@ARGV) != 2) {
    die "Usage: $0 vcfheader patient < tsv_from_suppl data S4"
}

my @H = ("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT");
my $vcfh = $ARGV[0];
my $patient = $ARGV[1];

sub start_vcf {
    my $name = shift;
    my $h = shift;
    my @this_h = @H;
    push(@this_h, $name);
    $name = 'trunk_' . $name . ".vcf";
    # put header!
    copy($h,$name) or die "Copy failed: $!";
    open my $fhandle,'>>', $name or die $!;
    print { $fhandle } join "\t", @this_h;
    print { $fhandle } "\n";
    return  $fhandle;
}

my $header = <STDIN>;
chomp $header;
my @fields = split("\t", $header);

my @keep_i = ();
my @vcf_files = ();
my $i = 0;
my %mrca = (
	P1 => ['17'],
	P2 => ['17_not_timed_to_WGD', '17_postWGD', '17_preWGD'],
	P3 => ['17']
);

# for trunk we have a single vcf for each patient
$vcf_files[0] = &start_vcf($patient, $vcfh);

foreach my $f (@fields) {
    if ($f =~ /($patient\.)[TN]\d\.\d_VAF/) {
        push(@keep_i, $i);
    }
    $i++;
}

sub find {
    my $wanted = shift;
    my $arref = shift;
    foreach my $af ( @{ $arref } ) {
        if ($af eq $wanted) {
            return 1;
        }
    }
    return 0;
}

while (<STDIN>) {
    chomp;
    my @line = split("\t", $_);
    my $chr = "chr".$line[0];
    my $pos = $line[1];
    my $ref = $line[2];
    my $alt = $line[3];
    my $id = $line[5];
    my $mrca = $line[4];
    print $mrca . "\n";
    my $qfif = "$chr\t$pos\t$id\t$ref\t$alt\t.\tPASS\tCONTQ=42\tGT:AF\t0/1:";
    my $tot = 0;
    for (my $k = 0; $k < scalar(@keep_i); $k++) { 
        # for trunk we print all muts found in at least 1 clone considered to be of the tumor MRCA
        if ($line[$keep_i[$k]] != 0 && &find($mrca, $mrca{$patient}) == 1) {
            $tot++;
        }
    }
    # we get the AF of the first clone
    if ($tot != 0) {
        print { $vcf_files[0] } $qfif . $line[$keep_i[0]] . "\n";
    }
}


close $_ foreach values @vcf_files;
