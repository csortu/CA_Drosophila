###########################################################################################
#	Program name:	fetch_cah_genomes.pl
#	Input:		"path for Directory" where log files are created by "get_cds_info.pl"
#	Usage:		fetch_cah_genomes.pl   `pwd`  (assuming script & log files are in cwd)
#	Output:		Parses log  files; Extract the accession numbers for genomic seqs &
#                     	downloads them to the "GENOMES" folder.
#	Author:		Preethy Nair
###########################################################################################

#!/usr/bin/perl -w
use diagnostics;
use strict;
use Bio::DB::GenBank;
use Bio::SeqIO;
use File::Path qw(make_path);
my $gb = new Bio::DB::GenBank;
chomp( @ARGV = <STDIN> ) unless @ARGV;
my $dir_for_log_files = $ARGV[0];
my @log_files         = <$dir_for_log_files/*.log>;
my $dir               = $dir_for_log_files . "/GENOMES/";

unless ( -d $dir ) {
    mkdir $dir or die "Can't create $dir:$!\n";
}

foreach my $log (@log_files) {
    open( my $f, '<', $log ) or die "Could not open file '$log' $!";
    while ( my $line = <$f> ) {
        chomp $line;
        my $nuc_acc = ( split( "\t", $line ) )[3];
        my $nuc_gb = $dir . $nuc_acc . ".gb";
        if ( !-e $nuc_gb ) {
            eval {
                my $seq    = $gb->get_Seq_by_acc($nuc_acc);
                my $seqout = Bio::SeqIO->new(
                    '-file'   => ">$nuc_gb",
                    '-format' => 'genbank'
                );
                $seqout->write_seq($seq);
            };
            print "$@ : $nuc_acc not found\n" if $@;
        }
    }
    close $f;
}

