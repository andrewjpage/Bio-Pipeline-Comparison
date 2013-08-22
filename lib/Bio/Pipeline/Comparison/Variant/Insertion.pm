package Bio::Pipeline::Comparison::Variant::Insertion;

# ABSTRACT: Represents a deletion

=head1 SYNOPSIS

Represents a deletion

    use Bio::Pipeline::Comparison::Variant::Insertion;
    my $obj = Bio::Pipeline::Comparison::Variant::Insertion->new( start_coord => 1, end_coord => 40, target_start_coord => 300);

=cut

use Moose;
use Bio::SeqIO;
use Bio::Pipeline::Comparison::Generate::VCFWriter;

has 'start_coord'           => ( is => 'ro', isa => 'Int', required => 1 );
has 'alternative_bases'     => ( is => 'ro', isa => 'Str', required => 1 );

no Moose;
__PACKAGE__->meta->make_immutable;
1;
