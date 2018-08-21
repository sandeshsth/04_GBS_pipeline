#!/usr/bin/perl
use warnings;
use strict;
use v5.16;

my $filename = $ARGV[0];

open IN, $filename or die "Cannot find file";
open OUT4, ">Count_Geno-$filename.txt";

my @count_2;
my $no_columns;
my @each_depth;
while (<IN>) {

	my $line = $_;

	# skip vcf header and keep heading with sample info
	if ( $_ =~ /^#/ ) {

		# put header in txt output
		if ( $_ =~ /^#CHROM/ ) {
			chomp $_;
			my @header = split( "\t", $_ );

			# number of columns
			$no_columns = scalar @header;
			for ( my $m = 9 ; $m < $no_columns ; $m++ ) {

			}
			print OUT4
"CHROM\tPOS\tREF\tALT\thet_count\tref_count\talt_count\tnull_count\ttotal_check\n";
		}
	}
	# finding DP and AD
	else {
		chomp $_;
		my @columns = split( "\t", $_ );

		# access format column to find DP and AD
		my @format_column = split( ":", $columns[8] );

		# If genotype has . then skip
		if ( $columns[3] eq "." or $columns[4] eq "." ) {
			next;
		}

		# If more than one allele
		if ( length $columns[3] > 1 or length $columns[4] > 1 ) {
			next;
		}

		# finding DP position:
		my $dp_position;
		my $allele_depth;

		# find (total depth in a sample) DP and AD positions
		for ( my $j = 0 ; $j < @format_column ; $j++ ) {
			if ( $format_column[$j] eq "DP" ) {
				$dp_position = $j;
			}
			if ( $format_column[$j] eq "AD" ) {
				$allele_depth = $j;
			}
		}

		my $call;
		my $ref_count  = 0;
		my $alt_count  = 0;
		my $het_count  = 0;
		my $null_count = 0;

		for ( my $i = 9 ; $i < $no_columns ; $i++ ) {
			my @each_column = split( ":", $columns[$i] );
			if (   $each_column[0] eq "0/0"
				or $each_column[0] eq "0/1"
				or $each_column[0] eq "1/0"
				or $each_column[0] eq "1/1" )
			{
				# Make the call
				# Hom ref
				if ( $each_column[0] eq "0/0" ) {
					$call = "$columns[3]$columns[3]";
					$ref_count++;
				}

				# Het
				elsif ( $each_column[0] eq "0/1" or $each_column[0] eq "1/0" ) {
					$call = "$columns[3]$columns[4]";
					$het_count++;
				}

				# Hom Alt
				elsif ( $each_column[0] eq "1/1" ) {
					$call = "$columns[4]$columns[4]";
					$alt_count++;
				}
			}

			# No genotype call in the VCF file
			else {
				$null_count++;
			}
		}
		my $total_check = $ref_count + $alt_count + $het_count + $null_count;
		print OUT4
"$columns[0]\t$columns[1]\t$columns[3]\t$columns[4]\t$het_count\t$ref_count\t$alt_count\t$null_count\t$total_check\n";
	}
}
say "Job Done, Please check txt file.";
exit;

