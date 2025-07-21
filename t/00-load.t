#-*-Perl-*-
use Test2::V0;
use IPC::Cmd qw(can_run);

my @required = qw(AnalyseDists RNAalifold RNAz R-scape);
my @missing  = grep { !can_run($_) } @required;

bail_out('Missing required external tools: ' . join(', ', @missing))
    if @missing;

use ok 'RNA';
diag( "Testing Vienna RNA $RNA::VERSION, Perl $], $^X" );

done_testing;
