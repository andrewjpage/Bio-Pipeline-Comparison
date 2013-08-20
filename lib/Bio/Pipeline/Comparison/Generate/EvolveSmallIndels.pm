
package Bio::Pipeline::Comparison::Generate::EvolveSmallIndels;

# ABSTRACT: Parse a fasta file and create small indels

=head1 SYNOPSIS

Parse a fasta file and create small indels
   use Bio::Pipeline::Comparison::Generate::EvolveSmallIndels;
   
   my $obj = Bio::Pipeline::Comparison::Generate::EvolveSmallIndels->new(
     fasta_file             => 'abc.fa',
     insertion_probability  => 0.0005,
     deletion_probability   => 0.0005
   );
   $obj->ids_to_gene_name;

=cut

use Moose;
use Bio::Tools::GFF;
use Bio::Pipeline::Comparison::Variant::Deletion;
use Bio::Pipeline::Comparison::Variant::Insertion;
use Math::Random qw( random_poisson );

has 'fasta_file'              => ( is => 'ro', isa => 'Str', required => 1 );
has '_fasta_parser'           => ( is => 'ro', isa => 'Bio::SeqIO', lazy    => 1, builder => '_build__fasta_parser' );
has 'insertion_probability'   => ( is => 'ro', isa => 'Int', default => 0.0005 );
has 'deletion_probability'    => ( is => 'ro', isa => 'Int', default => 0.0005 );
has 'indels'                  => ( is => 'ro', isa => 'ArrayRef',  lazy  => 1, builder => '_build_indels' );

sub _build_indels {
    my ($self) = @_;
    my @indels;

    while ( my $seq_obj = $self->_fasta_parser->next_seq() ) {
       my $sequence = $seq_obj->seq;
       
       for(my $i = 1; $i <= length($sequence); $i++ )
       {
        
        if($self->_has_there_been_a_random_hit($self->insertion_probability) == 1)
        {
          my $indel_length  = random_poisson(1, 3);
          
          
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

sub _build__fasta_parser {
    my ($self) = @_;
    return = Bio::SeqIO->new( -file => $self->fasta_file,         -format => 'Fasta' );
}


no Moose;
__PACKAGE__->meta->make_immutable;

1;
