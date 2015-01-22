########################################################################################################################
#	Program Name:	fetch_cah_prots.pl
# 	Usage:		fetch_cah_prots.pl  /INPUT_DIR/*
#			INPUT => Directory with (only) ortholog id files  (cah1_orts_id.txt .. ). 
#			Files provided by "Csaba Ortutay"
#
#	What it does:	Parses text files of CAH orthologs; picks Drosophila CAH orthologs
#			Gets "prot_acc" & "prot_gi" for each CAH ortholog from the file;
#                       Downnloads sequence file from GenBank; write them to respective CAH folders.
#       Author: 	Preethy Nair
#######################################################################################################################

#!/usr/bin/perl -w
use diagnostics;
use strict;
use Bio::DB::GenBank;
use Bio::SeqIO;
use File::Basename;
my $gb             = new Bio::DB::GenBank;
my @ortholog_files = @ARGV;

foreach my $file (@ortholog_files) {
    my $path = dirname($file);
    open( FILE, "$file" ) || print "Couldn't open file\n";
    my ( @prot_gi, @prot_acc, @species_name );
    my $dir = ( split( "_", basename($file) ) )[0];
    unless ( -d $dir ) {
        mkdir $dir or die "Can't create $dir:$!\n";
    }
    while (<FILE>) {
        chomp;
        if (/Drosophila/) {
            ( my $name, my $acc, my $gi ) = ( split("\t") )[ 0 .. 2 ];

            #   print "$name,$acc, $gi\n";
            push( @species_name, $name );
            push( @prot_acc,     $acc );
            push( @prot_gi,      $gi );
        }
    }

    #  download_prot_accs( \@prot_acc, $dir );
    download_prot_gis( \@prot_gi, $dir );
    close FILE;
}

sub download_prot_gis {

    #get sequences using gi & save
    my @prot_gis  = @{ $_[0] };
    my $directory = $_[1];
    foreach my $prot (@prot_gis) {
        eval {
            my $seq    = $gb->get_Seq_by_gi($prot);
            my $seqout = Bio::SeqIO->new(
                '-file'   => ">$directory/$prot.gb",
                '-format' => 'genbank'
            );
            $seqout->write_seq($seq);
        };
        print "$@ : $prot not found\n" if $@;
    }

}

sub download_prot_accs {

    #get sequences using acc & save
    my @protIds   = @{ $_[0] };
    my $directory = $_[1];
    print "@protIds\n";
    foreach my $prot (@protIds) {
        eval {
            my $seq    = $gb->get_Seq_by_acc($prot);
            my $seqout = Bio::SeqIO->new(
                '-file'   => ">$directory/$prot.gb",
                '-format' => 'genbank'
            );
            $seqout->write_seq($seq);
        };
        print "$@ : $prot not found \n" if $@;

    }
}
