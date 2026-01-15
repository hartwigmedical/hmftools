package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record SmallVariant(
        @NotNull DriverFields driver,
        @Nullable DriverCategory driverLikelihoodType,
        double driverLikelihood,
        @NotNull TranscriptImpact transcriptImpact,
        @Nullable TranscriptImpact otherImpact,
        boolean isCanonical,
        @NotNull PurpleVariantType type,
        @NotNull String gene,
        @NotNull String chromosome,
        int position,
        @NotNull String ref,
        @NotNull String alt,
        @NotNull PurpleCodingEffect worstCodingEffect,
        @NotNull HotspotType hotspot,
        @NotNull AllelicDepth tumorDepth,
        @Nullable AllelicDepth rnaDepth,
        double adjustedCopyNumber,
        double adjustedVAF,
        double minorAlleleCopyNumber,
        double variantCopyNumber,
        boolean biallelic,
        double biallelicProbability,
        @NotNull PurpleGenotypeStatus genotypeStatus,
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
            @NotNull Set<PurpleVariantEffect> effects,
            @NotNull PurpleCodingEffect codingEffect,
            boolean reported
    ) {}

    @RecordBuilder
    public record AllelicDepth(
            int totalReadCount,
            int alleleReadCount
    ) {}

    @NotNull @Override public String findingKey() { return driver.findingKey(); }
    @NotNull @Override public DriverSource driverSource() { return driver.driverSource(); }
    @NotNull @Override public ReportedStatus reportedStatus() { return driver.reportedStatus(); }
    @NotNull @Override public DriverInterpretation driverInterpretation() { return driver.driverInterpretation(); }

    public double clonalLikelihood()
    {
        return 1 - subclonalLikelihood;
    }
}
