package com.hartwig.hmftools.redux.splice;

// minimal lifted-coord summary for one primary alignment, cached in pass 1 so pass 2 can patch the
// partner record's mate fields (RNEXT / PNEXT / mate-strand / mate-unmapped / TLEN) without re-resolving.
// TODO: negativeStrand may need XS:A:+/- from gene strand for IGV colouring — revisit on redux integration.
public record LiftedMateInfo(
        String chromosome,
        int alignmentStart,
        int alignmentEnd,
        String liftedCigar,
        boolean negativeStrand,
        boolean unmapped)
{
    public static final LiftedMateInfo UNMAPPED = new LiftedMateInfo(null, 0, 0, null, false, true);

    public static LiftedMateInfo mapped(
            final String chromosome, final int alignmentStart, final int alignmentEnd, final String liftedCigar,
            final boolean negativeStrand)
    {
        return new LiftedMateInfo(chromosome, alignmentStart, alignmentEnd, liftedCigar, negativeStrand, false);
    }
}
