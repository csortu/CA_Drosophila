###############################################################################################
#	Script name:	make_cds.pl
#	Input:		(Directory name with) log files created by "get_cds_info.pl"
#	Usage:		make_cds.pl `pwd`  (when both the script & log files are in cwd)
#	Functions: 	picks accession for genomic files from log files;
#			makes cDNA sequence as per the annotation &
#               	saves it in fasta format to the CAH subfolders
#	Author:		Preethy Nair
################################################################################################
#!/usr/bin/perl -w
use diagnostics;
use strict;
use Bio::SeqIO;
use Bio::Tools::RestrictionEnzyme;
use File::Basename;
chomp( @ARGV = <STDIN> ) unless @ARGV;
my $dir_for_log_files = $ARGV[0];

#my @log_files         = <$dir_for_log_files/*.log>;

while ( defined( my $log = glob "${dir_for_log_files}/*.log" ) ) {
    print "$log\n";
    open( my $FH, '<', $log ) or die "Could not open file '$log' $!";
    while ( my $line = <$FH> ) {
        chomp($line);
        my $file_name = basename($log);
        my $ortho_dir = $file_name;
        $ortho_dir =~ s/\.log//;
        $ortho_dir =~ s/\_//;
        my $dirname = dirname($log);
        my ( $prot_acc, $dis_id, $location, $name, $nucleotide_acc );
        $prot_acc = ( split( "\t", $line ) )[2];
        $dis_id   = ( split( "\t", $line ) )[1];

        #Location contains the location row from the respective log files
        $location       = ( split( "\t", $line ) )[4];
        $name           = ( split( "\t", $line ) )[0];
        $nucleotide_acc = ( split( "\t", $line ) )[3];
        print "$nucleotide_acc\n";
        my $nuc_gb = "${dirname}/GENOMES/${nucleotide_acc}.gb";

        if ( -e $nuc_gb ) {
            my $cds = get_cds( $line, $nuc_gb, $location );

            #writing fasta files to the respective folders
            my $cds_seq = Bio::Seq->new(
                -seq        => $cds,
                -display_id => $dis_id,
                -desc       => "$name\t$prot_acc\t$nucleotide_acc\t$location",
                -alphabet   => "dna"
            );

            my $filename = "${dirname}/${ortho_dir}/${prot_acc}.fasta";
            my $seqout_obj =
              Bio::SeqIO->new( -file => ">$filename", -format => 'fasta' );
            $seqout_obj->write_seq($cds_seq);
        }
        else {
            print "$nuc_gb was not retieved \n";
        }
    }
    close $FH;
}

sub get_cds {

    my ( $gen_file_line, $nuc_file, $loc ) = @_;
    print "$loc\n";
    my $seqio_object =
      Bio::SeqIO->new( '-format' => 'GenBank', '-file' => "$nuc_file" );
    my $seq_obj    = $seqio_object->next_seq;
    my $id         = $seq_obj->display_id();
    my $seq_length = $seq_obj->length();
    my $sequence   = "";
    $_ = $loc;

    if ( $gen_file_line =~ /complement/ )
    {    # get exons + concatenate + reverse complement; if "complement" str.
        my $so;
        while (/(\d+)\s*\.\.\s*(\d+)/) {
            $loc = $';    #rest
            eval {
                my $string = $seq_obj->subseq( $1, $2 );
                $sequence = $sequence . $string;
                $so       = Bio::Seq->new(
                    -seq      => $sequence,
                    -alphabet => "dna"
                );
            };
            print "problem with $nuc_file $@" if $@;
            $_ = $loc;
        }

        $sequence = $so->primary_seq->revcom->seq;
    }
    else {    # get exons & concatenate; if not from complement strand
        while (/(\d+)\.\.(\d+)/) {
            eval {
                my $string = $seq_obj->subseq( $1, $2 );
                $sequence = $sequence . $string;
                $_        = $';
            };
            print "problem with $nuc_file $@" if $@;
        }
    }

    #the complete cds
    return ($sequence);
}
