#!/usr/bin/perl

use warnings;
use diagnostics;
use strict;

my @files = <*.cc>;
for (@files){
    if( /(qsim)(.*)/ ){
	#print "$2\n";
	my $new_name = 'MolPol' . $2;
	#print "$_ $new_name\n";
	rename $_, $new_name;
    }
}
