#########################################################################################
# Program Name:		get_cds_info.pl
# Usage:  		get_cds_info.pl  cah{1..9} cah{11..16}
#			cah{1..9} cah{11..16} are folders (with prot files) created by "fetch_cah_prots.pl"
# Output: 		Log file for each cah orthologs:
#			Information in each Log file:
#			"Species name\tprotein_accession\tprotein_gi\tnucleotide_accession\tcoded_by_row"
# Authors:		Preethy Nair and Csaba Ortutay
###################################################################################

#!/usr/bin/perl -w
use diagnostics;
use strict;
use Bio::SeqIO;
use Cwd 'abs_path';
my @cah_folders = @ARGV;

for my $dir (@cah_folders) {
    my $dir_path = abs_path($dir);
    my @gb_files = <$dir_path/*.gb>;
    my $log_file = $dir_path . ".log";
    open( my $FH, '>', $log_file ) or die "Could not open file '$log_file' $!";

    foreach my $gb_file (@gb_files) {
        my $seqio_object =
          Bio::SeqIO->new( '-format' => 'genbank', '-file' => $gb_file );
        my $seq_object     = $seqio_object->next_seq;
        my $prot_acc       = $seq_object->accession_number();
        my $primary_id     = $seq_object->primary_id();
        my $species_object = $seq_object->species;
        my $bi             = $species_object->binomial();
        my $coded_by_row;

        for my $feat_object ( $seq_object->get_SeqFeatures ) {
            if ( $feat_object->primary_tag eq "CDS" ) {
                if ( $feat_object->has_tag('coded_by') ) {
                    for my $val_1 ( $feat_object->get_tag_values('coded_by') ) {
                        $coded_by_row = $val_1;
                    }
                }
            }
        }

        # get DBSOURCE from $coded_by_row
        my $id = '';
        if ( $coded_by_row =~ /\(/ ) {
            $coded_by_row =~ /^.+\(([A-Z].+?):/;
            $id = $1;
        }
        else {
            $coded_by_row =~ /^([A-Z].+?):/;
            $id = $1;
        }
        print $FH "$bi\t $prot_acc\t$primary_id\t$id\t$coded_by_row\n";
    }
    close($FH);

    #  system("sort $log_file -o temp.out") == 0 or die ":$?";
    #  system("uniq temp.out >$log_file") == 0   or die ":$?";
    #  unlink "temp.out";
}
exit;

