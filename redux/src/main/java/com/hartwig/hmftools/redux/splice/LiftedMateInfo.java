package com.hartwig.hmftools.redux.splice;

// minimal lifted-coord summary for one primary alignment, cached in pass 1 so pass 2 can patch the
// partner record's mate fields (RNEXT / PNEXT / mate-strand / mate-unmapped / TLEN) without re-resolving.
public class LiftedMateInfo
{
    public final String Chromosome;
    public final int AlignmentStart;
    public final int AlignmentEnd;
    public final String LiftedCigar;
    public final boolean NegativeStrand; // TODO: IGV colouring: likely needs XS:A:+/- from gene strand, not this field; revisit
    public final boolean Unmapped;

    private LiftedMateInfo(
            final String chromosome, final int alignmentStart, final int alignmentEnd, final String liftedCigar,
            final boolean negativeStrand, final boolean unmapped)
    {
        Chromosome = chromosome;
        AlignmentStart = alignmentStart;
        AlignmentEnd = alignmentEnd;
        LiftedCigar = liftedCigar;
        NegativeStrand = negativeStrand;
        Unmapped = unmapped;
    }

    public static LiftedMateInfo mapped(
            final String chromosome, final int alignmentStart, final int alignmentEnd, final String liftedCigar,
            final boolean negativeStrand)
    {
        return new LiftedMateInfo(chromosome, alignmentStart, alignmentEnd, liftedCigar, negativeStrand, false);
    }

    public static LiftedMateInfo unmapped()
    {
        return new LiftedMateInfo(null, 0, 0, null, false, true);
    }
}
