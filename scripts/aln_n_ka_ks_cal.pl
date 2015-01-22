######################################################################################################
#	Script name:	aln_n_ka_ks_cal.pl
#	Usage:  	aln_n_ka_ks_cal.pl cah{1..9} cah{11..16}
#	       		cah ==> folders for each CAH ortholog with gb & fasta cDNA files
#       Function: 	Global pairwise protein alignment; DNA alignment;
#			dn, ds and Ka/Ks calculation & statistics
#	Output files:  	STAT_1.txt, ....., STAT_16.txt
#	Author:		Preethy Nair	
#######################################################################################################

#!/usr/bin/perl -w
use strict;
use diagnostics;
use Bio::SeqIO;
use Bio::Factory::EMBOSS;
use Bio::AlignIO;
use lib '$HOME/src/bioperl-live/Bio/Align/';
use Bio::Align::DNAStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Cwd;
use Cwd 'abs_path';
use File::Path qw(mkpath);

chomp( @ARGV = <STDIN> ) unless @ARGV;
my @cah_folders = @ARGV;

for my $dir (@cah_folders) {

    #print "$dir\n";

    my $dir_path = abs_path($dir);

    # cds (fasta) and protein (gb )files
    my @file_fasta   = <$dir_path/*.fasta>;
    my @file_genbank = <$dir_path/*.gb>;
    my $cwd          = getcwd;

    #dir for the statistics
    my $stat_dir = "${cwd}/STAT";
    unless ( -d $stat_dir ) {

        eval { mkpath($stat_dir) };
        print "Couldn't create $stat_dir: $@" if $@;

    }
    my $no = $dir;    #ortholog no:
    $no =~ s/cah//;
    my $stat_file = "${stat_dir}/STAT_${no}.txt";

    open( my $FH, ">", $stat_file );
    print $FH "SPECIES_NO:1\tProt_acc_1\tSPECIES_NO:2\t",
      "Prot_acc_2\tD_n\tD_s\tZ-SCORE\tKa_Ks\n";
    my $n    = @file_fasta;    #The total number of fasta files
    my $loop = 0;

    # loop through each pair of orthologs
    for ( my $i = 0 ; $i < $n ; $i++ ) {
        for ( my $j = $i + 1 ; $j < $n ; $j++ ) {
            $loop++;

            #protein no:1 & its cds file
            my $n_file_1 = $file_fasta[$i];
            my $prot_file_1 = ( split( /\./, $n_file_1 ) )[0] . ".gb";

            #protein no:2 & its cds file
            my $n_file_2 = $file_fasta[$j];
            my $prot_file_2 = ( split( /\./, $n_file_2 ) )[0] . ".gb";

            #protein pairwise-global alignment & other details
            my $prot_info = aln_aa( $prot_file_1, $prot_file_2, $no );

            # parsing the dna files
            my ( $seqio_dna_1, $seqio_dna_2, $seq_dna_1, $seq_dna_2 );
            $seqio_dna_1 =
              Bio::SeqIO->new( '-format' => 'Fasta', '-file' => $n_file_1 );
            $seqio_dna_2 =
              Bio::SeqIO->new( '-format' => 'Fasta', '-file' => $n_file_2 );
            $seq_dna_1 = $seqio_dna_1->next_seq;
            $seq_dna_2 = $seqio_dna_2->next_seq;

            # checking for nucleotide translation
            my ( $dna_transl_1, $prot_1, $n_id_1, $dna_transl_2, $prot_2,
                $n_id_2 );

            $dna_transl_1 = ( $seq_dna_1->translate() )->seq();
            $prot_1       = $prot_info->{prot_seq_1};
            chop($dna_transl_1);
            $n_id_1 = $seq_dna_1->display_id();

            $dna_transl_2 = ( $seq_dna_2->translate() )->seq();
            $prot_2       = $prot_info->{prot_seq_2};
            chop($dna_transl_2);
            $n_id_2 = $seq_dna_2->display_id();

            if (   ( $dna_transl_1 eq $prot_1 )
                && ( $dna_transl_2 eq $prot_2 ) )
            {

                # dna hash
                my %dna = (
                    $prot_info->{pr_id_1} => $seq_dna_1,
                    $prot_info->{pr_id_2} => $seq_dna_2
                );
                my $dna_aln_dir = "${cwd}/DnaAlign";
                unless ( -d $dna_aln_dir ) {

                    mkpath($dna_aln_dir)
                      or die "Can't create $dna_aln_dir:$!\n";
                }

                # aminoacid alignment to nucleotide space
                my $aa_aln = $prot_info->{aln_obj};
                my $dna_aln = &aa_to_dna_aln( $aa_aln, \%dna );
                my ( $seq_dna_id, $seq_dna2_id ) =
                  map { $_->display_id } $dna_aln->each_seq;

                my $dna_alignout = new Bio::AlignIO(
                    -format => 'msf',
                    -file   => ">>${dna_aln_dir}/dna_aln_CAH_${no}.msf"
                );

                $dna_alignout->write_aln($dna_aln);

                # creating the statistics
                my $stats = new Bio::Align::DNAStatistics;
                my $results =
                  $stats->calc_KaKs_pair( $dna_aln, $seq_dna_id, $seq_dna2_id );
                my @k  = keys %{ $results->[0] };
                my $nk = @k;

                #print "key count: $nk\n";
                if ( $nk > 0 ) {
                    my $seq_1 = $results->[0]{Seq1};
                    my $seq_2 = $results->[0]{Seq2};
                    my $dn    = $results->[0]{D_n};
                    my $ds    = $results->[0]{D_s};
                    my $z     = $results->[0]{z_score};
                    my $y     = $dn / $ds;

                    print $FH "$prot_info->{binomial_1}\t",
                      "$prot_info->{pr_acc_1}\t",
                      "$prot_info->{binomial_2} \t",
                      " $prot_info->{pr_acc_2}\t$dn\t$ds\t$z\t$y\n";
                }
                else {

                    #  print "Problem Pairs: $n_file_1, $n_file_2\n";
                    #exit;
                }
            }
            else {

                # print "Protein & Nucleotide translations not same\n";
                # print "\n$n_id_1\n$n_id_2\n";
                # exit;
            }

        }
    }
    close($FH);
}

sub aln_aa {

    my ( $prt_file_1, $prt_file_2, $cah_no ) = @_;
    my ( $seqio_prot_1, $seqio_prot_2, $seq_prot_1, $seq_prot_2, $prot_acc_1,
        $prot_acc_2, $bi_1, $bi_2 );
    $seqio_prot_1 =
      Bio::SeqIO->new( '-format' => 'genbank', '-file' => $prt_file_1 );
    $seqio_prot_2 =
      Bio::SeqIO->new( '-format' => 'genbank', '-file' => $prt_file_2 );

    $seq_prot_1 = $seqio_prot_1->next_seq;
    $seq_prot_2 = $seqio_prot_2->next_seq;
    $prot_acc_1 = $seq_prot_1->accession_number();
    $prot_acc_2 = $seq_prot_2->accession_number();
    $bi_1       = $seq_prot_1->species->binomial();
    $bi_2       = $seq_prot_2->species->binomial();

    my $seq1 = $seq_prot_1->seq();
    my $seq2 = $seq_prot_2->seq();

    my $sp_1 = $bi_1;
    $sp_1 =~ s/^(\w+\s)//;
    $sp_1 = substr( $sp_1, 0, 6 );
    my $sp_2 = $bi_2;
    $sp_2 =~ s/^(\w+\s)//;
    $sp_2 = substr( $sp_2, 0, 6 );
    my $current_dir = getcwd;
    my $needle_dir  = "${current_dir}/Needle";

    #dir for the statistics
    my $aa_aln_dir = "${current_dir}/ProtAlign";
    unless ( -d $aa_aln_dir ) {

        mkpath($aa_aln_dir) or die "Can't create $aa_aln_dir:$!\n";
    }

    unless ( -d $needle_dir ) {
        mkpath($needle_dir) or die "Can't create $needle_dir:$!\n";
    }

    # aligning the protein sequences using emboss needle
    my ( $factory, $needle, $needle_out, $alnin, $aln, $alignout );
    $factory = new Bio::Factory::EMBOSS;
    $needle  = $factory->program("needle");

    $needle_out = "${needle_dir}/${sp_1}_${sp_2}.needle";

    $needle->run(
        {
            -asequence => $seq_prot_1,
            -bsequence => $seq_prot_2,
            -outfile   => $needle_out
        }
    );
    $alnin    = new Bio::AlignIO( -format => 'emboss', -file => $needle_out );
    $aln      = $alnin->next_aln;
    $alignout = new Bio::AlignIO(
        -format => 'msf',
        -file   => ">${aa_aln_dir}/${sp_1}_${sp_2}.msf"
    );
    $alignout->write_aln($aln);

    #all aln
    my $alignout_all = new Bio::AlignIO(
        -format => 'msf',
        -file   => ">>${aa_aln_dir}/Prot_Aln_Cah_${cah_no}.msf"
    );
    $alignout_all->write_aln($aln);
    my ( $seqid, $seq2id ) = map { $_->display_id } $aln->each_seq;
    my $prot_info_ref = {
        pr_id_1    => $seqid,
        pr_id_2    => $seq2id,
        prot_seq_1 => $seq1,
        prot_seq_2 => $seq2,
        pr_acc_1   => $prot_acc_1,
        pr_acc_2   => $prot_acc_2,
        aln_obj    => $aln,
        binomial_1 => $bi_1,
        binomial_2 => $bi_2
    };
    return ($prot_info_ref);
}
