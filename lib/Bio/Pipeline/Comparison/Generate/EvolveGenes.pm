
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

=cut

use Moose;
use Bio::Tools::GFF;
use Bio::Pipeline::Comparison::Variant::Deletion;
use Bio::Pipeline::Comparison::Variant::Insertion;

has 'gff_file'              => ( is => 'ro', isa => 'Str',             required => 1 );
has '_gff_parser'           => ( is => 'ro', isa => 'Bio::Tools::GFF', lazy     => 1, builder => '_build__gff_parser' );
has 'insertion_probability' => ( is => 'ro', isa => 'Num',             default  => 0.0005 );
has 'deletion_probability'  => ( is => 'ro', isa => 'Num',             default  => 0.0005 );
has 'inversion_probability' => ( is => 'ro', isa => 'Num',             default  => 0.0005 );

has 'indels'     => ( is => 'ro', isa => 'ArrayRef', lazy => 1, builder => '_build_indels' );
has '_all_cdses' => ( is => 'ro', isa => 'ArrayRef', lazy => 1, builder => '_build__all_cdses' );
has '_randomise_cdses_index' =>
  ( is => 'ro', isa => 'ArrayRef', lazy => 1, builder => '_build__randomise_cdses_index' );

has '_total_length' => ( is => 'rw', isa => 'Int', default => 0 );

sub _build__all_cdses {
    my ($self) = @_;
    my @cdses;
    my $total_length = 0;
    while ( my $segment = $self->_gff_parser->next_segment ) {
        my $segment_length = $segment->end;
        $total_length += $segment_length;
        while ( my $raw_feature = $self->_gff_parser->next_feature() ) {
            last unless defined($raw_feature);    # No more features
            next if !( $raw_feature->primary_tag eq 'CDS' );
            push( @cdses, $raw_feature );
        }
    }
    $self->_total_length($total_length);
    return \@cdses;
}

sub _build__randomise_cdses_index {
    my ($self)         = @_;
    my @all_cdses      = @{ $self->_all_cdses };
    my @shuffled_index = shuffle( 0 .. $#all_cdses );
    return \@shuffled_index;
}

sub _delete_genes {
    my ( $self, $indels ) = @_;
    my $num_genes_to_delete = int( $self->deletion_probability * $self->_total_length );
    for ( my $i = 0 ; $i < $num_genes_to_delete ; $i++ ) {
        my $cdsindex = pop( @{ $self->_randomise_cdses_index } );
        push(
            @{$indels},
            Bio::Pipeline::Comparison::Variant::Deletion->new(
                start_coord => $self->_all_cdses->[$cdsindex]->start,
                end_coord   => $self->_all_cdses->[$cdsindex]->end
            )
        );
    }
}

sub _insert_genes {
    my ( $self, $indels ) = @_;
    my $num_genes_to_insert = int( $self->insertion_probability * $self->_total_length );
    for ( my $i = 0 ; $i < $num_genes_to_insert ; $i++ ) {
        my $cdsindex           = pop( @{ $self->_randomise_cdses_index } );
        my $size_of_gene       = $self->_all_cdses->[$cdsindex]->end - $self->_all_cdses->[$cdsindex]->start;
        my $target_start_coord = ( int( rand( $self->_total_length - $size_of_gene ) ) );

        if ( $target_start_coord == 0 || ( $target_start_coord + $size_of_gene > $self->_total_length ) ) {
            $target_start_coord = 1;
        }

        push(
            @{$indels},
            Bio::Pipeline::Comparison::Variant::GeneInsertion->new(
                start_coord        => $self->_all_cdses->[$cdsindex]->start,
                end_coord          => $self->_all_cdses->[$cdsindex]->end,
                target_start_coord => $target_start_coord
            )
        );
    }
}

sub _invert_genes {
    my ( $self, $indels ) = @_;
    my $num_genes_to_invert = int( $self->inversion_probability * $self->_total_length );
    for ( my $i = 0 ; $i < $num_genes_to_invert ; $i++ ) {
        my $cdsindex = pop( @{ $self->_randomise_cdses_index } );
        push(
            @{$indels},
            Bio::Pipeline::Comparison::Variant::Inversion->new(
                start_coord => $self->_all_cdses->[$cdsindex]->start,
                end_coord   => $self->_all_cdses->[$cdsindex]->end
            )
        );
    }
}

sub _build_indels {
    my ($self) = @_;
    my @indels;

    $self->_delete_genes( \@indels );
    $self->_invert_genes( \@indels );
    $self->_insert_genes( \@indels );

    return \@indels;
}

sub _build__gff_parser {
    my ($self) = @_;
    open( my $fh, $self->gff_file ) or die "Couldnt open GFF file";
    my $gff_parser = Bio::Tools::GFF->new( -fh => $fh, gff_version => 3 );
    return $gff_parser;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
