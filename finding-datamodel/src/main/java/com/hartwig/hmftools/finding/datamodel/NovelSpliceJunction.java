package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record NovelSpliceJunction(
        @NotNull String findingKey,
        @NotNull String gene,
        @NotNull String chromosome,
        int junctionStart,
        int junctionEnd,
        @NotNull Type type,
        int exonStart,
        int exonEnd,
        int fragmentCount,
        int depthStart,
        int depthEnd,
        @NotNull Context regionStart,
        @NotNull Context regionEnd,
        int cohortFrequency
) implements Finding
{
    public enum Type
    {
        SKIPPED_EXONS,
        NOVEL_5_PRIME,
        NOVEL_3_PRIME,
        NOVEL_EXON,
        NOVEL_INTRON,
        MIXED_TRANS,
        INTRONIC,
        EXON_INTRON,
        CIRCULAR,
        UNKNOWN
    }

    public enum Context
    {
        SPLICE_JUNC,
        EXONIC,
        INTRONIC,
        MIXED,
        UNKNOWN
    }
}
