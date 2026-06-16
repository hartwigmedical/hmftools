package com.hartwig.hmftools.tars.liftback;

public enum LiftBackCategory
{
    UNMAPPED,
    SUPPLEMENTARY,

    // position falls in an inter-transcript N-spacer or otherwise outside any transcript span
    LIFT_FAILED,

    // ---- 1 locus ----

    REF_SINGLE,
    TX_SINGLE,

    // same locus + CIGAR, no N - both views agree
    BOTH_AGREE,

    // tx has N-junction, ref is soft-clipped - tx wins (picks up splice ref missed)
    BOTH_TX_JUNCTION_REF_SOFTCLIP,

    // tx has N-junction, ref is a full match through the intron with low NM - ref wins (unspliced/retained intron)
    BOTH_TX_JUNCTION_REF_MATCH,

    // tx soft-clipped at exon boundary, ref full match - ref wins (likely intron retention)
    BOTH_TX_SOFTCLIP_REF_MATCH,

    // single locus, no discriminator rule fires
    BOTH_AMBIGUOUS,

    // ---- >=2 loci ----

    REF_MULTI,
    TX_MULTI,

    // distinct loci; tx has N-junction, ref alts don't (intronless paralogs) - tx is biologically faithful
    BOTH_MULTI_TX_JUNCTION,

    // distinct loci, tx-junction rule didn't fire - genuine multi-mapper (paralog/family overlap)
    BOTH_MULTI;

    public PrimaryBucket primaryBucket()
    {
        return switch(this)
        {
            case REF_SINGLE, TX_SINGLE, BOTH_AGREE, BOTH_TX_JUNCTION_REF_SOFTCLIP, BOTH_TX_JUNCTION_REF_MATCH,
                    BOTH_TX_SOFTCLIP_REF_MATCH, BOTH_MULTI_TX_JUNCTION -> PrimaryBucket.RESOLVED;
            case TX_MULTI, REF_MULTI, BOTH_MULTI, BOTH_AMBIGUOUS -> PrimaryBucket.PROBLEM;
            default -> PrimaryBucket.NON_PRIMARY;
        };
    }

    public enum PrimaryBucket
    {
        RESOLVED,
        PROBLEM,
        NON_PRIMARY
    }
}
