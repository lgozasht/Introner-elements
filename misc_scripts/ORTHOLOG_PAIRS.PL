#!/usr/bin/perl
($#ARGV == 1) || die "Usage: ORTHOLOG_PAIRS.PL sp1.v.sp2 sp2.v.sp1";







foreach $file (@ARGV) {

    open (IN, $file);

    

    while (<IN>) {

	$already{(split)[0]}++ ||

	    $hits{join ("\t", sort ((split)[0,1]))}++;

    }

}



foreach $key (keys %hits) {

    ($hits{$key} == 2) && 

	print "$key\n";

}