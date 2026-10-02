package com.hartwig.hmftools.tars.liftback;

public record AlignmentKey(String chromosome, int position, String cigar, boolean forwardStrand)
{
    public static AlignmentKey from(final LiftedAlignment alignment)
    {
        return new AlignmentKey(
                alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar, alignment.ForwardStrand);
    }

    public Locus locus()
    {
        return new Locus(chromosome, position);
    }

    public record Locus(String chromosome, int position)
    {
    }
}
