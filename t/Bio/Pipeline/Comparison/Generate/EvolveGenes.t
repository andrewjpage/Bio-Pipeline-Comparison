#!/usr/bin/env perl
use strict;
use warnings;
BEGIN { unshift( @INC, './lib' ) }

BEGIN {
    use Test::Most;
    use_ok('Bio::Pipeline::Comparison::Generate::EvolveGenes');
    use_ok('Bio::Pipeline::Comparison::Variant::Deletion');
    use_ok('Bio::Pipeline::Comparison::Variant::Insertion');
    use_ok('Bio::Pipeline::Comparison::Variant::GeneInsertion');
    use_ok('Bio::Pipeline::Comparison::Variant::Inversion');
}

ok(my $obj = Bio::Pipeline::Comparison::Generate::EvolveGenes->new(gff_file => 't/data/example_annotation.gff',), 'initialise object');


done_testing();