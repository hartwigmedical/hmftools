package com.hartwig.hmftools.redux.splice;

// per-record decision categories assigned by LiftBackResolver.
public enum LiftBackCategory
{
    // input record was unmapped on entry; passed through.
    UNMAPPED,

    // supplementary records take a separate bucket; not run through the discriminator tree.
    SUPPLEMENTARY,

    // record's primary alignment was on a tx contig but ContigTranslator failed (e.g. position falls in
    // an inter-transcript N-spacer of the packed alt contig [case to be made for supplementary]). Emitted unmapped with placeholder coords.
    LIFT_FAILED,

    // ---- 1 locus ----

    // record's lifted alignment set contains only ref-chromosome alignments collapsing to a single locus.
    REF_SINGLE,

    // tx-contig alignments only (no ref alts), all collapsing to a single genomic locus after lift.
    TX_SINGLE,

    // ref + tx both contributed, single locus, same CIGAR, no N operator. Confirms the read's alignment
    // is identical whether viewed through the genome or a transcript.
    BOTH_AGREE,

    // single locus with both ref and tx alignments. Tx has an N (junction) in its CIGAR; ref is soft-clipped
    // around the junction. The transcript-aware alignment wins, picks up the splice the genome alignment missed.
    BOTH_TX_JUNCTION_REF_SOFTCLIP,

    // single locus with both ref and tx. Tx has an N (junction); ref alignment is a full match across the
    // supposed intron with no soft-clip. Ref full-match with low NM is overwhelming evidence the read is
    // genuinely unspliced (pre-mRNA / retained-intron / DNA contamination) — the tx N-CIGAR is the
    // artifact, not the ref full-match. We pick ref and drop the tx alt.
    BOTH_TX_JUNCTION_REF_MATCH,

    // single locus with both ref and tx alignments. Tx has soft-clip at an exon boundary; ref is a full match.
    // Ref alignment wins — likely an intron-retention event the tx-contig encoding can't capture.
    BOTH_TX_SOFTCLIP_REF_MATCH,

    // single locus with both ref and tx but none of the discriminator rules above fire. Falls through.
    BOTH_AMBIGUOUS,

    // ---- >=2 loci ----

    // ref-chromosome alignments only (no tx alts) spanning >=2 distinct genomic loci.
    // Genuine multi-mapper on the reference alone, paralog or repeat hit.
    REF_MULTI,

    // tx-contig alignments only but lift maps them to >=2 distinct genomic loci. Multi-locus, ambiguous.
    TX_MULTI,

    // ref + tx both contributed at distinct loci. Tx has an N (real annotated junction); none of the ref
    // alts have an N (they sit at intronless paralogs — typically processed pseudogenes). The tx alignment
    // is the biologically faithful one. We promote tx to primary and demote ref alts.
    BOTH_MULTI_TX_JUNCTION,

    // ref + tx both contributed but they map to >=2 distinct genomic loci, and the tx-junction discriminator
    // didn't fire. Genuine multi-mapper — paralog hit, gene-family overlap, or similar.
    BOTH_MULTI;

    // headline grouping for at-a-glance stats / TSV scanning. Most reads fall in RESOLVED on a typical sample;
    // PROBLEM is the bucket worth digging into (multi-locus + discriminator outcomes); NON_PRIMARY is everything
    // the resolver didn't make a positive decision about.
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
                // discriminator outcomes are confident picks — bucket as RESOLVED, not PROBLEM
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

    public enum PrimaryBucket
    {
        // a single locus identified with confidence — the read is downstream-ready.
        RESOLVED,
        // multi-locus or a discriminator outcome — downstream tools need to be aware which sub-category this came from.
        PROBLEM,
        // not a primary alignment the resolver made a decision about (unmapped, lift failed, or supplementary).
        NON_PRIMARY
    }
}
