package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record RnaFusion(
        @NotNull String findingKey,
        @Nullable String geneStart,
        @Nullable String geneEnd,
        @NotNull String chromosomeStart,
        @NotNull String chromosomeEnd,
        int positionStart,
        int positionEnd,
        @NotNull String junctionTypeStart,
        @NotNull String junctionTypeEnd,
        @NotNull KnownType knownType,
        @NotNull StructuralVariantType structuralVariantType,
        int splitFragments,
        int realignedFragments,
        int discordantFragments,
        int depthStart,
        int depthEnd,
        int cohortFrequency
) implements Finding
{
    public enum KnownType
    {
        NONE,
        PROMISCUOUS_3,
        PROMISCUOUS_5,
        PROMISCUOUS_BOTH,
        KNOWN_PAIR,
        EXON_DEL_DUP
    }

    public enum StructuralVariantType
    {
        BND,
        DEL,
        DUP,
        INF,
        INS,
        INV,
        SGL
    }

    @NotNull
    public String display()
    {
        return String.format("%s::%s", geneStart != null ? geneStart : "", geneEnd != null ? geneEnd : "");
    }
}
