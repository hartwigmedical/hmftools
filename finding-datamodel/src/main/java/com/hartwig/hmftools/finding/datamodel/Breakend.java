package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record Breakend(
    int id,
    int svId,
    @NotNull GeneOrientation geneOrientation,
    boolean disruptive,
    @NotNull TranscriptRegionType regionType,
    @NotNull TranscriptCodingType codingType,
    int nextSpliceExonRank,
    int orientation,
    int exonUp,
    int exonDown,
    double junctionCopyNumber)
{
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
