package Bio::Pipeline::Comparison::Variant::Deletion;

# ABSTRACT: Represents a deletion

=head1 SYNOPSIS

Represents a deletion

    use Bio::Pipeline::Comparison::Variant::Deletion;
    my $obj = Bio::Pipeline::Comparison::Variant::Deletion->new(start_coord => 1, end_coord => 20);

=cut

use Moose;
use Bio::SeqIO;
use Bio::Pipeline::Comparison::Generate::VCFWriter;

has 'start_coord'  => ( is => 'ro', isa => 'Int', required => 1 );
has 'end_coord'    => ( is => 'ro', isa => 'Int', required => 1 );

no Moose;
__PACKAGE__->meta->make_immutable;
1;
