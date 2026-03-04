package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record Breakend(
    int id,
    int svId,
    @NotNull String gene,
    @NotNull String chromosome,
    @NotNull String chromosomeBand,
    @NotNull String transcript,
    boolean isCanonical,
    @NotNull GeneOrientation geneOrientation,
    boolean disruptive,
    boolean reported,
    double undisruptedCopyNumber,
    @NotNull Type type,
    @NotNull TranscriptRegionType regionType,
    @NotNull TranscriptCodingType codingType,
    int nextSpliceExonRank,
    int orientation,
    int exonUp,
    int exonDown,
    double junctionCopyNumber)
{
    public enum Type
    {
        BND,
        DEL,
        DUP,
        INF,
        INS,
        INV,
        SGL
    }

    public enum GeneOrientation
    {
        UPSTREAM,
        DOWNSTREAM
    }

    public enum TranscriptRegionType
    {
        UNKNOWN,
        UPSTREAM,
        EXONIC,
        INTRONIC,
        IG,
        DOWNSTREAM
    }

    public enum TranscriptCodingType
    {
        UNKNOWN,
        CODING,
        UTR_5P,
        UTR_3P,
        NON_CODING,
        ENHANCER
    }
}
