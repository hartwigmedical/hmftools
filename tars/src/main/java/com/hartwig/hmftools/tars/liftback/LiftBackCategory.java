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

    // same locus + CIGAR, no N — both views agree
    BOTH_AGREE,

    // tx has N-junction, ref is soft-clipped — tx wins (picks up splice ref missed)
    BOTH_TX_JUNCTION_REF_SOFTCLIP,

    // tx has N-junction, ref is a full match through the intron with low NM — ref wins (unspliced/retained intron)
    BOTH_TX_JUNCTION_REF_MATCH,

    // tx soft-clipped at exon boundary, ref full match — ref wins (likely intron retention)
    BOTH_TX_SOFTCLIP_REF_MATCH,

    // single locus, no discriminator rule fires
    BOTH_AMBIGUOUS,

    // ---- >=2 loci ----

    REF_MULTI,
    TX_MULTI,

    // distinct loci; tx has N-junction, ref alts don't (intronless paralogs) — tx is biologically faithful
    BOTH_MULTI_TX_JUNCTION,

    // distinct loci, tx-junction rule didn't fire — genuine multi-mapper (paralog/family overlap)
    BOTH_MULTI;

    public PrimaryBucket primaryBucket()
    {
        switch(this)
        {
            case REF_SINGLE:
            case TX_SINGLE:
            case BOTH_AGREE:
            case BOTH_TX_JUNCTION_REF_SOFTCLIP:
            case BOTH_TX_JUNCTION_REF_MATCH:
            case BOTH_TX_SOFTCLIP_REF_MATCH:
            case BOTH_MULTI_TX_JUNCTION:
                return PrimaryBucket.RESOLVED;
            case TX_MULTI:
            case REF_MULTI:
            case BOTH_MULTI:
            case BOTH_AMBIGUOUS:
                return PrimaryBucket.PROBLEM;
            case UNMAPPED:
            case LIFT_FAILED:
            case SUPPLEMENTARY:
            default:
                return PrimaryBucket.NON_PRIMARY;
        }
    }

    // tx-contig won over bwa; lift-back has strictly better placement info so the STAR ladder runs unconditionally.
    public boolean txWon()
    {
        switch(this)
        {
            case TX_SINGLE:
            case TX_MULTI:
            case BOTH_TX_JUNCTION_REF_SOFTCLIP:
            case BOTH_MULTI_TX_JUNCTION:
                return true;
            default:
                return false;
        }
    }

    public enum PrimaryBucket
    {
        RESOLVED,
        PROBLEM,
        NON_PRIMARY
    }
}
