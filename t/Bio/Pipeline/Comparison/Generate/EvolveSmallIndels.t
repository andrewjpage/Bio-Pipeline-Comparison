#!/usr/bin/env perl
use strict;
use warnings;
BEGIN { unshift( @INC, './lib' ) }

BEGIN {
    use Test::Most;
    use_ok('Bio::Pipeline::Comparison::Generate::EvolveSmallIndels');
    use_ok('Bio::Pipeline::Comparison::Variant::Deletion');
    use_ok('Bio::Pipeline::Comparison::Variant::Insertion');
}


ok(my $obj = Bio::Pipeline::Comparison::Generate::EvolveSmallIndels->new(
   fasta_file             => 't/data/test_genome.fa',
   insertion_probability  => 0.01,
   deletion_probability   => 0.01
 ), 'initalise a test genome');

ok(my $all_indel_objs = $obj->indels, 'extract out all the variant objects');

ok(@{$all_indel_objs} > 10 && @{$all_indel_objs} < 30, 'check the right number are created' );

my $count_insertion = 0;
my $count_deletion = 0;

for my $indel_obj (@{$all_indel_objs})
{
  ok(($indel_obj->isa('Bio::Pipeline::Comparison::Variant::Deletion') || $indel_obj->isa('Bio::Pipeline::Comparison::Variant::Insertion') ), 'Check the right types are coming out');
  if($indel_obj->isa('Bio::Pipeline::Comparison::Variant::Deletion'))
  {
    $count_deletion++;
    ok($indel_obj->start_coord < $indel_obj->end_coord, 'start coord must be less than end coord');
  }

  if($indel_obj->isa('Bio::Pipeline::Comparison::Variant::Insertion'))
  {
    $count_insertion++;
    ok(length($indel_obj->alternative_bases() ) > 0, 'alternative bases must be present');
    ok($indel_obj->alternative_bases() =~ /^[ACGTacgt]+$/, 'Only standard bases allowed');
  }
}

ok($count_deletion > 0 && $count_deletion < 20, 'number of deletions is correct');
ok($count_insertion > 0 && $count_insertion < 20, 'number of insertions is correct');

my $del1 = Bio::Pipeline::Comparison::Variant::Deletion->new( start_coord => 1, end_coord => 4);
my $insertion1 = Bio::Pipeline::Comparison::Variant::Insertion->new( start_coord => 6, alternative_bases => 'ACGGG');
ok(my $obj_overlapping = Bio::Pipeline::Comparison::Generate::EvolveSmallIndels->new(
   fasta_file             => 't/data/test_genome.fa',
   insertion_probability  => 0.01,
   deletion_probability   => 0.01,
   indels => [$del1, $insertion1]
 ), 'initalise a test genome with overlapping indels');
is( $obj_overlapping->_is_overlapping(Bio::Pipeline::Comparison::Variant::Insertion->new( start_coord => 1, alternative_bases => 'A')), 1, 'check overlap at start of deletion' );
is( $obj_overlapping->_is_overlapping(Bio::Pipeline::Comparison::Variant::Insertion->new( start_coord => 2, alternative_bases => 'A')), 1, 'check overlap in middle of deletion' );
is( $obj_overlapping->_is_overlapping(Bio::Pipeline::Comparison::Variant::Insertion->new( start_coord => 4, alternative_bases => 'A')), 0, 'check nonoverlap at end of deletion' );
is( $obj_overlapping->_is_overlapping(Bio::Pipeline::Comparison::Variant::Deletion->new( start_coord => 5, end_coord => 7)), 1, 'check overlap of deletion with insertion' );
is( $obj_overlapping->_is_overlapping(Bio::Pipeline::Comparison::Variant::Deletion->new( start_coord => 5, end_coord => 6)), 0, 'check nonoverlap of deletion beside insertion' );

done_testing();