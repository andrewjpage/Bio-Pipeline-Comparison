
package Bio::Pipeline::Comparison::Generate::EvolveGenes;

# ABSTRACT: 

=head1 SYNOPSIS

Parse a GFF and identify target genes for insertions and deletions
   use Bio::Pipeline::Comparison::Generate::EvolveGenes;
   
   my $obj = Bio::Pipeline::Comparison::Generate::EvolveGenes->new(
     gff_file               => 'abc.gff',
     insertion_probability  => 0.0005,
     deletion_probability   => 0.0005
   );
   $obj->ids_to_gene_name;

=cut

use Moose;
use Bio::Tools::GFF;
use Bio::Pipeline::Comparison::Variant::Deletion;
use Bio::Pipeline::Comparison::Variant::Insertion;

has 'gff_file'              => ( is => 'ro', isa => 'Str', required => 1 );
has '_gff_parser'           => ( is => 'ro', isa => 'Bio::Tools::GFF', lazy    => 1, builder => '_build__gff_parser' );
has 'insertion_probability' => ( is => 'ro', isa => 'Int', default => 0.0005 );
has 'deletion_probability'  => ( is => 'ro', isa => 'Int', default => 0.0005 );

has 'indels'                => ( is => 'ro', isa => 'ArrayRef',  lazy  => 1, builder => '_build_indels' );

sub _build_indels {
    my ($self) = @_;
    my @indels;

    while (my $segment = $self->_gff_parser->next_segment)
    {
      my $segment_length = $segment->end;
      while ( my $raw_feature = $self->_gff_parser->next_feature() ) {
          last unless defined($raw_feature);    # No more features
          next if !( $raw_feature->primary_tag eq 'CDS' );
          my @junk;
          
          if($self->_has_there_been_a_random_hit($self->insertion_probability) == 1)
          {
            my $target_start_coord = (int(rand($segment_length))+ $raw_feature->end)% $segment_length;
            if ($target_start_coord == 0) 
            {
              $target_start_coord = 1;
            }
            push(@indels, Bio::Pipeline::Comparison::Variant::Insertion->new( start_coord => $raw_feature->start, end_coord => $raw_feature->end, target_start_coord => $target_start_coord));
          }
          elsif($self->_has_there_been_a_random_hit($self->deletion_probability) == 1)
          {
            push(@indels, Bio::Pipeline::Comparison::Variant::Deletion->new(start_coord => $raw_feature->start, end_coord => $raw_feature->end));
          }
      }
    }
    
    $self->_gff_parser->close();
    return \@indels;
}


sub _has_there_been_a_random_hit
{
  my ($self, $probability) = @_;
  if ( rand(1) <= $probability ) {
    return 1;
  }
  return 0;
}

sub _build__gff_parser {
    my ($self) = @_;
    open( my $fh,  $self->gff_file ) or die "Couldnt open GFF file";
    my $gff_parser = Bio::Tools::GFF->new( -fh => $fh, gff_version => 3 );
    return $gff_parser;
}


no Moose;
__PACKAGE__->meta->make_immutable;

1;
