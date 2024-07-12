package com.hartwig.hmftools.orange.conversion;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGermlineDeletion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineDetectionMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.orange.algo.purple.CodingEffectDeterminer;

import org.jetbrains.annotations.NotNull;

public final class PurpleConversion
{
    @NotNull
    public static PurpleCopyNumber convert(@NotNull com.hartwig.hmftools.common.purple.PurpleCopyNumber copyNumber)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(copyNumber.chromosome())
                .start(copyNumber.start())
                .end(copyNumber.end())
                .averageTumorCopyNumber(copyNumber.averageTumorCopyNumber())
                .build();
    }

    @NotNull
    public static PurpleGeneCopyNumber convert(@NotNull GeneCopyNumber geneCopyNumber)
    {
        return ImmutablePurpleGeneCopyNumber.builder()
                .gene(geneCopyNumber.geneName())
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .minMinorAlleleCopyNumber(geneCopyNumber.minMinorAlleleCopyNumber())
                .build();
    }

    @NotNull
    public static PurpleDriver convert(@NotNull DriverCatalog catalog)
    {
        return ImmutablePurpleDriver.builder()
                .gene(catalog.gene())
                .transcript(catalog.transcript())
                .type(PurpleDriverType.valueOf(catalog.driver().name()))
                .driverLikelihood(catalog.driverLikelihood())
                .likelihoodMethod(PurpleLikelihoodMethod.valueOf(catalog.likelihoodMethod().name()))
                .isCanonical(catalog.isCanonical())
                .build();
    }

    @NotNull
    public static PurpleQC convert(@NotNull com.hartwig.hmftools.common.purple.PurpleQC purpleQC)
    {
        return ImmutablePurpleQC.builder()
                .status(ConversionUtil.mapToIterable(purpleQC.status(), PurpleConversion::convert))
                .germlineAberrations(ConversionUtil.mapToIterable(purpleQC.germlineAberrations(), PurpleConversion::convert))
                .amberMeanDepth(purpleQC.amberMeanDepth())
                .contamination(purpleQC.contamination())
                .totalCopyNumberSegments(purpleQC.copyNumberSegments())
                .unsupportedCopyNumberSegments(purpleQC.unsupportedCopyNumberSegments())
                .deletedGenes(purpleQC.deletedGenes())
                .build();
    }

    @NotNull
    public static PurpleAllelicDepth convert(@NotNull AllelicDepth allelicDepth)
    {
        return ImmutablePurpleAllelicDepth.builder()
                .totalReadCount(allelicDepth.TotalReadCount)
                .alleleReadCount(allelicDepth.AlleleReadCount)
                .build();
    }

    @NotNull
    public static PurpleGermlineDeletion convert(@NotNull GermlineDeletion germlineDeletion)
    {
        return ImmutablePurpleGermlineDeletion.builder()
                .gene(germlineDeletion.GeneName)
                .chromosome(germlineDeletion.Chromosome)
                .chromosomeBand(germlineDeletion.ChromosomeBand)
                .regionStart(germlineDeletion.RegionStart)
                .regionEnd(germlineDeletion.RegionEnd)
                .depthWindowCount(germlineDeletion.DepthWindowCount)
                .exonStart(germlineDeletion.ExonStart)
                .exonEnd(germlineDeletion.ExonEnd)
                .detectionMethod(PurpleGermlineDetectionMethod.valueOf(germlineDeletion.DetectionMethod.name()))
                .normalStatus(PurpleGermlineStatus.valueOf(germlineDeletion.NormalStatus.name()))
                .tumorStatus(PurpleGermlineStatus.valueOf(germlineDeletion.TumorStatus.name()))
                .germlineCopyNumber(germlineDeletion.GermlineCopyNumber)
                .tumorCopyNumber(germlineDeletion.TumorCopyNumber)
                .filter(germlineDeletion.Filter)
                .cohortFrequency(germlineDeletion.CohortFrequency)
                .reported(germlineDeletion.Reported)
                .build();
    }

    @NotNull
    public static PurpleGermlineAberration convert(@NotNull GermlineAberration aberration)
    {
        return PurpleGermlineAberration.valueOf(aberration.name());
    }

    @NotNull
    public static PurpleQCStatus convert(@NotNull com.hartwig.hmftools.common.purple.PurpleQCStatus qcStatus)
    {
        return PurpleQCStatus.valueOf(qcStatus.name());
    }

    @NotNull
    public static PurpleCodingEffect convert(@NotNull CodingEffect effect)
    {
        return PurpleCodingEffect.valueOf(effect.name());
    }

    @NotNull
    public static PurpleVariantEffect convert(@NotNull VariantEffect effect)
    {
        return PurpleVariantEffect.valueOf(effect.name());
    }

    @NotNull
    public static PurpleTranscriptImpact convert(@NotNull VariantTranscriptImpact impact, boolean reported)
    {
        List<VariantEffect> effectsList = VariantEffect.effectsToList(impact.Effects);
        List<PurpleVariantEffect> purpleEffects = ConversionUtil.mapToList(effectsList, PurpleConversion::convert);
        PurpleCodingEffect purpleCodingEffect = convert(CodingEffectDeterminer.determineCodingEffect(effectsList));

        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(impact.Transcript)
                .hgvsCodingImpact(impact.HgvsCoding)
                .hgvsProteinImpact(impact.HgvsProtein)
                .inSpliceRegion(impact.SpliceRegion)
                .effects(purpleEffects)
                .codingEffect(purpleCodingEffect)
                .reported(reported)
                .build();
    }
}
