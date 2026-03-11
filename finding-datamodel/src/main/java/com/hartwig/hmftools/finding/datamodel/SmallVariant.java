package com.hartwig.hmftools.finding.datamodel;

import java.util.List;
import java.util.Set;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record SmallVariant(
        @NotNull DriverFields driver,
        @Nullable DriverCategory driverLikelihoodType,
        @NotNull TranscriptImpact transcriptImpact,
        @Nullable TranscriptImpact otherImpact,
        boolean isCanonical,
        @NotNull VariantType type,
        @NotNull String gene,
        @NotNull String chromosome,
        int position,
        @NotNull String ref,
        @NotNull String alt,
        @NotNull CodingEffect worstCodingEffect,
        @NotNull HotspotType hotspot,
        @NotNull AllelicDepth allelicDepth,
        @Nullable AllelicDepth rnaDepth,
        double adjustedCopyNumber,
        double adjustedVAF,
        double minorAlleleCopyNumber,
        double variantCopyNumber,
        boolean biallelic,
        double biallelicLikelihood,
        @NotNull GenotypeStatus genotypeStatus,
        int repeatCount,
        double subclonalLikelihood,
        @Nullable List<Integer> localPhaseSets
) implements Driver
{
    @RecordBuilder
    public record TranscriptImpact(
            @NotNull String transcript,
            @NotNull String hgvsCodingImpact,
            @NotNull String hgvsProteinImpact,
            @Nullable Integer affectedCodon,
            @Nullable Integer affectedExon,
            boolean inSpliceRegion,
            @NotNull Set<VariantEffect> effects,
            @NotNull CodingEffect codingEffect,
            boolean reported
    )
    {
    }

    @RecordBuilder
    public record AllelicDepth(
            int totalReadCount,
            int alleleReadCount
    )
    {
    }

    public enum VariantType
    {
        MNP,
        SNP,
        INDEL,
        UNDEFINED
    }

    public enum CodingEffect
    {
        NONSENSE_OR_FRAMESHIFT,
        SPLICE,
        MISSENSE,
        SYNONYMOUS,
        NONE,
        UNDEFINED
    }

    public enum VariantEffect
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
        OTHER
    }

    public enum HotspotType
    {
        HOTSPOT,
        NEAR_HOTSPOT,
        NON_HOTSPOT
    }

    public enum GenotypeStatus
    {
        HOM_REF,
        HET,
        HOM_ALT,
        UNKNOWN
    }

    @NotNull
    @Override
    public String findingKey()
    {
        return driver.findingKey();
    }

    @NotNull
    @Override
    public DriverSource driverSource()
    {
        return driver.driverSource();
    }

    @NotNull
    @Override
    public ReportedStatus reportedStatus()
    {
        return driver.reportedStatus();
    }

    @NotNull
    @Override
    public DriverInterpretation driverInterpretation()
    {
        return driver.driverInterpretation();
    }

    @Override
    public double driverLikelihood()
    {
        return driver.driverLikelihood();
    }

    @Override
    public boolean isReported()
    {
        return driver.isReported();
    }

    @Override
    public Set<String> genes()
    {
        return Set.of(gene());
    }

    public double clonalLikelihood()
    {
        return 1 - subclonalLikelihood;
    }
}
