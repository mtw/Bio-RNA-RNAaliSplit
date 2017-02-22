#!perl -T
use 5.010;
use strict;
use warnings;
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'Bio::RNA::AliSplit' ) || print "Bail out!\n";
}

diag( "Testing Bio::RNA::AliSplit $Bio::RNA::AliSplit::VERSION, Perl $], $^X" );
