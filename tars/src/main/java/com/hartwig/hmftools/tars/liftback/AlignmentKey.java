package com.hartwig.hmftools.tars.liftback;

public record AlignmentKey(String chromosome, int position, String cigar, boolean forwardStrand)
{
    public static AlignmentKey from(final LiftedAlignment alignment)
    {
        return new AlignmentKey(
                alignment.liftedPosition().Chromosome, alignment.liftedPosition().Position, alignment.liftedCigar(),
                alignment.forwardStrand());
    }

    public Locus locus()
    {
        return new Locus(chromosome, position);
    }

    public record Locus(String chromosome, int position)
    {
    }
}
