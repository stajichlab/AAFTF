#!/usr/bin/env perl

use File::Spec;
use strict;
use warnings;

my %stats;
my $model = 'capnodiales_odb10';
my $BUSCO_dir = 'BUSCO';
my $prefix = shift || 'Rant';
my $telomere_report = 'telomere_reports';
my $dir = shift || ".";
my @header;
my %header_seen;

opendir(DIR,$dir) || die $!;
my $first = 1;

for my $subdir ( readdir(DIR) ) {

	my $file = "$prefix.sorted.stats.txt";
	my $filepath = File::Spec->catfile($dir,$subdir,$file);
	if ( -d "$dir/$subdir" && -f $filepath ) {
		my $stem = join(".",$subdir,$prefix);

		open(my $fh => $filepath) || die "cannot open $filepath: $!";
		while(<$fh>) {
			next if /^\s+$/;
			s/^\s+//;
			chomp;
			if ( /\s*(.+)\s+=\s+(\d+(\.\d+)?)/ ) {
				my ($name,$val) = ($1,$2);
				$name =~ s/\s*$//;
				$name =~ s/\s+/_/g;
				$stats{$stem}->{$name} = $val;

				if( ! $header_seen{$name} ) {
					push @header, $name;
					$header_seen{$name} = 1;
				}
			}
		}

		if ( -d $telomere_report ) {
			if ( $first ) {
				push @header, qw(Telomeres_Found Telomeres_Fwd Telomeres_Rev Telomeres_CompleteChrom);
			}
			my $telomerefile = File::Spec->catfile($telomere_report,sprintf("%s.sorted.telomere_report.txt",$stem));

			if ( -f $telomerefile ) {
				open(my $fh => $telomerefile) || die $!;
				my %contigs_with_tel;
				while(<$fh>) {
					if( /^(\S+)\s+(forward|reverse)\s+(\S+)/i ){
						$contigs_with_tel{$1}->{$2} = $3;
					} elsif (/^Telomeres found:\s+(\d+)\s+\((\S+)\s+forward,\s+(\S+)\s+reverse\)/ ) {
						$stats{$stem}->{'Telomeres_Found'} = $1;
						$stats{$stem}->{'Telomeres_Fwd'} = $2;
						$stats{$stem}->{'Telomeres_Rev'} = $3;
					}
				}
				for my $ctg ( keys %contigs_with_tel ) {
					if (exists $contigs_with_tel{$ctg}->{'forward'} &&
						exists $contigs_with_tel{$ctg}->{'reverse'} ) {
						$stats{$stem}->{'Telomeres_CompleteChrom'} +=1; # or ++ but count up the number of times we have a ctg w fwd&rev
					}
				}
			}

		}

		if ( -d $BUSCO_dir ) {
			if ( $first ) {
				push @header, qw(BUSCO_Complete BUSCO_Single BUSCO_Duplicate
								 BUSCO_Fragmented BUSCO_Missing BUSCO_NumGenes
					);
			}

			my $busco_file = File::Spec->catfile($BUSCO_dir,$subdir,
												 sprintf("short_summary.specific.%s.%s.txt",$model,$subdir));


			if ( -f $busco_file ) {
				open(my $fh => $busco_file) || die $!;
				while(<$fh>) {
					if (/^\s+C:(\d+\.\d+)\%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)/ ) {
						$stats{$stem}->{"BUSCO_Complete"} = $1;
						$stats{$stem}->{"BUSCO_Single"} = $2;
						$stats{$stem}->{"BUSCO_Duplicate"} = $3;
						$stats{$stem}->{"BUSCO_Fragmented"} = $4;
						$stats{$stem}->{"BUSCO_Missing"} = $5;
						$stats{$stem}->{"BUSCO_NumGenes"} = $6;
					}
				}

			} else {
				warn("Cannot find $busco_file");
			}
		}
		$first = 0;
	}
}
print join("\t", qw(SampleID), @header), "\n";
foreach my $sp ( sort keys %stats ) {
	print join("\t", $sp, map { $stats{$sp}->{$_} || 'NA' } @header), "\n";
}
