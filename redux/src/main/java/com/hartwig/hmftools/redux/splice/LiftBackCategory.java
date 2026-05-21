package com.hartwig.hmftools.redux.splice;

public enum LiftBackCategory
{
    UNMAPPED,
    SUPPLEMENTARY,

    // primary alignment was on a tx contig but ContigTranslator couldn't lift it (e.g. position falls in
    // an inter-transcript N-spacer of the packed alt contig).
    LIFT_FAILED,

    // ---- 1 locus ----

    REF_SINGLE,
    TX_SINGLE,

    // ref + tx both contributed, same locus + CIGAR, no N. The two views of the read agree.
    BOTH_AGREE,

    // tx has an N (junction); ref is soft-clipped around it. Tx wins — picks up the splice ref missed.
    BOTH_TX_JUNCTION_REF_SOFTCLIP,

    // tx has an N (junction); ref is a full match through the supposed intron with low NM. Ref wins —
    // overwhelming evidence the read is genuinely unspliced (pre-mRNA / retained-intron / DNA contamination).
    BOTH_TX_JUNCTION_REF_MATCH,

    // tx soft-clipped at an exon boundary; ref is a full match. Ref wins — likely intron retention the
    // tx-contig encoding can't represent.
    BOTH_TX_SOFTCLIP_REF_MATCH,

    // both ref + tx at single locus, no discriminator rule fires.
    BOTH_AMBIGUOUS,

    // ---- >=2 loci ----

    REF_MULTI,
    TX_MULTI,

    // ref + tx at distinct loci. Tx has an N (annotated junction); ref alts have no N (intronless paralogs,
    // typically processed pseudogenes). Tx is the biologically faithful pick.
    BOTH_MULTI_TX_JUNCTION,

    // ref + tx at distinct loci, tx-junction rule didn't fire. Genuine multi-mapper — paralog / family overlap.
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

    // categories where the discriminator promoted a tx-contig alignment over what bwa would have picked.
    // In these cases the lift-back has strictly better placement info than bwa, so the STAR ladder runs
    // unconditionally; elsewhere it caps at the input MAPQ.
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
