#!/usr/bin/perl -w 
use strict;

my $file=$ARGV[0];

open(IN,$file);
my @lines=<IN>;
close IN;

chomp @lines;

my $ofile=$file;

$ofile=~s/\.txt/_ka_ks.txt/;

open (OUTF,">$ofile");

foreach my $li (@lines){

    my @parts=split(/\t/,$li);

    map{$_=~s/^\s*//}@parts;
    map{$_=~s/\s*$//}@parts;

#    print "$li\n";

    my $newli=join("\t",@parts);
    print OUTF "$newli\n";


}

close OUTF;


exit;

