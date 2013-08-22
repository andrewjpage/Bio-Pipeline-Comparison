
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
has 'insertion_probability'   => ( is => 'ro', isa => 'Num', default => 0.0005 );
has 'deletion_probability'    => ( is => 'ro', isa => 'Num', default => 0.0005 );
has 'indels'                  => ( is => 'ro', isa => 'ArrayRef',  lazy  => 1, builder => '_build_indels' );


sub _next_coord
{
  my ($self, $sequence_length, $current_coord) = @_;
  my $total_probab = $self->insertion_probability + $self->deletion_probability;
  my $number_of_variants =  $sequence_length * $total_probab;
  my $step_size = int($sequence_length / $number_of_variants);
  my $next_increment =  $step_size;
  if(rand(2)>= 1)
  {
    $next_increment  = $next_increment - int(rand(int($step_size*0.05)));
  }
  else
  {
    $next_increment  = $next_increment + int(rand(int($step_size*0.05)));
  }
  return $next_increment;
}

sub _build_indels {
    my ($self) = @_;
    my @indels;

    while ( my $seq_obj = $self->_fasta_parser->next_seq() ) {
       my $sequence = $seq_obj->seq;
       
       my $sequence_length = length($sequence);
       for(my $i = 1; $i <= $sequence_length; $i += $self->_next_coord($sequence_length) )
       {
        
          if(rand(2)>= 1)
          {
            my $indel_length  =  $self->_indel_length;
            
            my $target_start_coord = (int(rand($sequence_length))+ $i)% $sequence_length;
            if ($target_start_coord == 0 || ($target_start_coord + $indel_length >= $sequence_length)) 
            {
              $target_start_coord = 1;
            }
            
            my $proposed_variant =  Bio::Pipeline::Comparison::Variant::Insertion->new( start_coord => $i, alternative_bases => substr( $sequence, $target_start_coord-1, $indel_length));
            if($self->_is_overlapping($proposed_variant) == 0)
            {
              push(@indels, $proposed_variant );
            }
          }
          else{
        
            my $indel_length  = $self->_indel_length;
            my $proposed_variant = Bio::Pipeline::Comparison::Variant::Deletion->new(start_coord => $i, end_coord => ($i + $indel_length));
            if($self->_is_overlapping($proposed_variant) == 0)
            {
              push(@indels, $proposed_variant);
            }
        }
      }
    }
    
    $self->_fasta_parser->close();
    return \@indels;
}

sub _indel_length
{
   my ($self) = @_;
   my $indel_length  = int(random_poisson(1, 3));
   if($indel_length < 1)
   {
     return 1;
   }
   else
   {
     return $indel_length ;
   }
}

sub _is_overlapping
{
  my ($self, $variation_obj) = @_;
  
  if($variation_obj->isa('Bio::Pipeline::Comparison::Variant::Insertion'))
  {
    return 1 if($self->_is_coord_overlapping($variation_obj->start_coord) == 1);
  }
  elsif($variation_obj->isa('Bio::Pipeline::Comparison::Variant::Deletion'))
  {
    for( my $deletion_coord = $variation_obj->start_coord; $deletion_coord < $variation_obj->end_coord; $deletion_coord++)
    {
      return 1 if( $self->_is_coord_overlapping($deletion_coord) == 1);
    }
  }
  return 0;
}

sub _is_coord_overlapping
{
   my ($self, $coord) = @_;
   for my $indel (@{$self->indels})
   {
      if($indel->isa('Bio::Pipeline::Comparison::Variant::Insertion'))
      {
        return 1 if($indel->start_coord == $coord );
      }
      elsif($indel->isa('Bio::Pipeline::Comparison::Variant::Deletion'))
      {
        for( my $deletion_coord = $indel->start_coord; $deletion_coord < $indel->end_coord; $deletion_coord++)
        {
          return 1 if($deletion_coord == $coord);
        }
      }
   }
   return 0;
   
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
    return  Bio::SeqIO->new( -file => $self->fasta_file,         -format => 'Fasta' );
}


no Moose;
__PACKAGE__->meta->make_immutable;

1;
