package com.hartwig.hmftools.tars.liftback;

// Which contig source the final primary came from. The discriminator's decision, reduced to three
// values; the reason rides alongside in DecidingFeature and the locus count in LiftBackResult.numLoci.
public enum Outcome
{
    REF,
    TX,
    UNRESOLVED
}
