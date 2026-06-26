package com.hartwig.hmftools.tars.liftback;

// Why the primary was resolved the way it was - the evidence that decided ref vs tx. Null for the
// non-decision record states (unmapped / lift-failed / supplementary), which carry no contest.
public enum DecidingFeature
{
    SOLE_REF,                // only ref aligned
    SOLE_TX,                 // only tx aligned
    CONCORDANT,              // ref and tx agree, both contiguous, no junction
    REF_READS_THROUGH,       // ref reads contiguously where tx carries a junction - read is unspliced
    JUNCTION,                // tx spans a trusted junction the ref cannot - tx wins
    JUNCTION_OVER_CONTIGUOUS,  // multi-locus: tx junction beats intronless ref alts
    MULTIMAPPER,             // multi-locus, no junction to break the tie
    AMBIGUOUS;               // single locus, no rule fires

    // The outcome is a strict function of the feature, so it is derived here rather than stored alongside.
    public Outcome outcome()
    {
        return switch(this)
        {
            case SOLE_REF, CONCORDANT, REF_READS_THROUGH -> Outcome.REF;
            case SOLE_TX, JUNCTION, JUNCTION_OVER_CONTIGUOUS -> Outcome.TX;
            case MULTIMAPPER, AMBIGUOUS -> Outcome.UNRESOLVED;
        };
    }
}
