#!/usr/bin/perl -w
use strict;

foreach my $i (1..9,13..16,11..12){
    my $inf='/home/ortutay/Dokumentumok/Sync/Munka/KaKs/STAT_'.$i.'_ka_ks.txt';
    next if not -f $inf;

    my $gene;
    if ($i == 11) {
	$gene='CARP-A';
    }elsif($i == 12){
	$gene='CARP-B';
    }else{
	$gene='CAH'.$i;
    }

    print "Statistics for $gene orthologs\n";

    my @lines=qx/cat $inf/;

    print @lines;
    print "\n\n";

}

exit;

