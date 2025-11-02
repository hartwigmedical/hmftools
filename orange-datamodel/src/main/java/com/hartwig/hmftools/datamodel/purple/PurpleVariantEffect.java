package com.hartwig.hmftools.datamodel.purple;

public enum PurpleVariantEffect
{
    STOP_GAINED,
    STOP_LOST,
    START_LOST,
    FRAMESHIFT,
    SPLICE_ACCEPTOR,
    SPLICE_DONOR,
    INFRAME_INSERTION,
    INFRAME_DELETION,
    MISSENSE,
    PHASED_MISSENSE,
    PHASED_INFRAME_INSERTION,
    PHASED_INFRAME_DELETION,
    SYNONYMOUS,
    PHASED_SYNONYMOUS,
    INTRONIC,
    FIVE_PRIME_UTR,
    THREE_PRIME_UTR,
    UPSTREAM_GENE,
    NON_CODING_TRANSCRIPT,
    OTHER;

    public String effect()
    {
        return switch (this) {
            case STOP_GAINED -> "stop_gained";
            case STOP_LOST -> "stop_lost";
            case START_LOST -> "start_lost";
            case FRAMESHIFT -> "frameshift_variant";
            case SPLICE_ACCEPTOR -> "splice_acceptor_variant";
            case SPLICE_DONOR -> "splice_donor_variant";
            case INFRAME_INSERTION -> "inframe_insertion";
            case INFRAME_DELETION -> "inframe_deletion";
            case MISSENSE -> "missense_variant";
            case PHASED_INFRAME_INSERTION -> "phased_inframe_insertion";
            case PHASED_INFRAME_DELETION -> "phased_inframe_deletion";
            case PHASED_MISSENSE -> "phased_missense";
            case SYNONYMOUS -> "synonymous_variant";
            case PHASED_SYNONYMOUS -> "phased_synonymous";
            case INTRONIC -> "intron_variant";
            case FIVE_PRIME_UTR -> "5_prime_UTR_variant";
            case THREE_PRIME_UTR -> "3_prime_UTR_variant";
            case UPSTREAM_GENE -> "upstream_gene_variant";
            case NON_CODING_TRANSCRIPT -> "non_coding_transcript_exon_variant";
            default -> "other";
        };
    }
}
