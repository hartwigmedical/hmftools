package com.hartwig.hmftools.orange.conversion;

import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
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
import com.hartwig.hmftools.datamodel.purple.PurpleSomaticLikelihood;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.orange.algo.purple.CodingEffectDeterminer;

import org.jetbrains.annotations.NotNull;

public final class PurpleConversion
{
    public static PurpleCopyNumber convert(final com.hartwig.hmftools.common.purple.PurpleCopyNumber copyNumber)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(copyNumber.chromosome())
                .start(copyNumber.start())
                .end(copyNumber.end())
                .averageTumorCopyNumber(copyNumber.averageTumorCopyNumber())
                .build();
    }

    public static PurpleGeneCopyNumber convert(final GeneCopyNumber geneCopyNumber)
    {
        return ImmutablePurpleGeneCopyNumber.builder()
                .gene(geneCopyNumber.geneName())
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.ChromosomeBand)
                .transcript(geneCopyNumber.TransName)
                .isCanonical(geneCopyNumber.IsCanonical)
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .maxCopyNumber(geneCopyNumber.maxCopyNumber())
                .minMinorAlleleCopyNumber(geneCopyNumber.MinMinorAlleleCopyNumber)
                .build();
    }

    public static PurpleDriver convert(final DriverCatalog catalog)
    {
        return ImmutablePurpleDriver.builder()
                .gene(catalog.gene())
                .transcript(catalog.transcript())
                .type(PurpleDriverType.valueOf(catalog.driver().name()))
                .driverLikelihood(catalog.driverLikelihood())
                .driverInterpretation(DriverInterpretation.interpret(catalog.driverLikelihood()))
                .likelihoodMethod(PurpleLikelihoodMethod.valueOf(catalog.likelihoodMethod().name()))
                .isCanonical(catalog.isCanonical())
                .reportedStatus(com.hartwig.hmftools.datamodel.driver.ReportedStatus.valueOf(catalog.reportedStatus().name()))
                .build();
    }

    public static PurpleQC convert(final com.hartwig.hmftools.common.purple.PurpleQC purpleQC)
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

    public static PurpleAllelicDepth convert(final AllelicDepth allelicDepth)
    {
        return ImmutablePurpleAllelicDepth.builder()
                .totalReadCount(allelicDepth.TotalReadCount)
                .alleleReadCount(allelicDepth.AlleleReadCount)
                .build();
    }

    public static PurpleGermlineDeletion convert(final GermlineAmpDel germlineDeletion)
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
                .normalStatus(convert(germlineDeletion.NormalStatus))
                .tumorStatus(convert(germlineDeletion.TumorStatus))
                .germlineCopyNumber(germlineDeletion.GermlineCopyNumber)
                .tumorCopyNumber(germlineDeletion.TumorCopyNumber)
                .filter(germlineDeletion.Filter)
                .cohortFrequency(germlineDeletion.CohortFrequency)
                .reported(germlineDeletion.Reported == ReportedStatus.REPORTED)
                .build();
    }

    public static PurpleGermlineAberration convert(final GermlineAberration aberration)
    {
        return PurpleGermlineAberration.valueOf(aberration.name());
    }

    public static PurpleQCStatus convert(final com.hartwig.hmftools.common.purple.PurpleQCStatus qcStatus)
    {
        return PurpleQCStatus.valueOf(qcStatus.name());
    }

    public static PurpleCodingEffect convert(final CodingEffect effect)
    {
        return PurpleCodingEffect.valueOf(effect.name());
    }

    public static PurpleVariantEffect convert(final VariantEffect effect)
    {
        return PurpleVariantEffect.valueOf(effect.name());
    }

    public static PurpleTranscriptImpact convert(final VariantTranscriptImpact impact, boolean reported)
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

    public static PurpleGermlineStatus convert(final GermlineStatus germlineStatus)
    {
        return PurpleGermlineStatus.valueOf(germlineStatus.name());
    }

    public static PurpleSomaticLikelihood convert(final com.hartwig.hmftools.common.variant.SomaticLikelihood somaticLikelihood)
    {
        return PurpleSomaticLikelihood.valueOf(somaticLikelihood.name());
    }
}
