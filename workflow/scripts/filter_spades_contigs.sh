perl - ${snakemake_params[opts]} > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}" <<"EOF"
#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

my $usage =
"Usage: $0 --csv slimmed_assembly_graph.csv --paths contigs.paths --fasta contigs.fasta\n";
my ( $graph_csv, $contigs_path, $contigs_fasta, $help );

GetOptions(
    'csv=s'   => \$graph_csv,
    'paths=s' => \$contigs_path,
    'fasta=s' => \$contigs_fasta,
    'help|h'  => \$help
) or die $usage;

if ($help) {
    print $usage;
    exit 0;
}

print {*STDERR}
"INPUT:\n\tCSV: $graph_csv\n\tPaths: $contigs_path\n\tFASTA: $contigs_fasta\n";

# collect all edges belonging to a target
open( my $fh_csv, '<', $graph_csv )
  or die "Could not open file '$graph_csv' $!";
my %edges;
while (<$fh_csv>) {
    my ( $edge, $target ) = split /\t|,/, $_;
    next if $edge eq 'EDGE';
    $edges{$edge} = $target;
}
close $fh_csv;

print {*STDERR} "edges: ", join( ", ", sort { $a <=> $b } keys %edges ), "\n";

# Read contig paths and collect edges in each contig
open( my $fh_contigs, '<', $contigs_path )
  or die "Could not open file '$contigs_path' $!";

# my %paths;
my %contigs;

# my @edge_list;
my $id;
while (<$fh_contigs>) {
    s/\R//g;    # Remove newline characters
    if (/^NODE_.*/) {
        $id = $_;
    }
    else {
        my @split = map { s/(-|\+);?//; $_ } split /,/, $_;
        for my $e (@split) {
            push @{ $contigs{$id} }, $e;
        }
    }
}
close $fh_contigs;

# Check which contigs are completely represented by the organelle regions
my ( @matched, @unmatched, @partial );
my %contig_types;
my %filter;
for my $contig ( sort sort_node_ids keys %contigs ) {
    my @edge_list = @{ $contigs{$contig} };
    my $match     = 0;
    my %types;
    for my $e (@edge_list) {
        if ( $edges{$e} ) {
            $match++;
            $types{ $edges{$e} }++;
        }
    }
    $contig_types{$contig} = \%types;
    if ( $match == scalar(@edge_list) ) {
        push @matched, $contig;
        $filter{$contig} = 1;
    }
    elsif ( $match == 0 ) {
        push @unmatched, $contig;
    }
    else {
        push @partial, $contig;
    }
}

# Log contigs with matches
print {*STDERR} "Matched contigs:\n";
for my $contig ( sort sort_node_ids @matched ) {
    next if $contig =~ /'$/;
    my @types = keys %{ $contig_types{$contig} };
    print {*STDERR} "\t$contig\t", join( ", ", @types ), "\n";
}
print {*STDERR} "\n";

# print {*STDERR} "Unmatched contigs:\n", join("\n\t", @unmatched), "\n\n";
print {*STDERR} "Partially matched contigs:\n";
for my $contig ( sort sort_node_ids @partial ) {
    next if $contig =~ /'$/;
    my @types = keys %{ $contig_types{$contig} };
    print {*STDERR} "\t$contig\t", join( ", ", @types ), "\n";
}
print {*STDERR} "\n";

# Filter the contigs fasta file based on the matched contigs
open( my $fh_fasta, '<', $contigs_fasta )
  or die "Could not open file '$contigs_fasta' $!";
my $print;
while (<$fh_fasta>) {
    if (/^>(.*)/) {
        my $id = $1;
        if ( $filter{$id} ) {
            $print = 0;
        }
        else {
            $print = 1;
        }
    }
    print $_ if $print;
}
close $fh_fasta;

sub sort_node_ids {
    my ($a_id) = $a =~ /NODE_(\d+)_/;
    my ($b_id) = $b =~ /NODE_(\d+)_/;
    return $a_id <=> $b_id;
}
EOF
