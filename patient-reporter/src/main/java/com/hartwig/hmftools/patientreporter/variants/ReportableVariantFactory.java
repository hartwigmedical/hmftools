package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.germline.DriverGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.somatic.DriverSomaticVariant;

import org.jetbrains.annotations.NotNull;

final class ReportableVariantFactory {

    private ReportableVariantFactory() {
    }

    @NotNull
    static List<ReportableVariant> mergeSomaticAndGermlineVariants(@NotNull List<DriverSomaticVariant> somaticVariantsReport,
            @NotNull List<DriverGermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingChoice) {
        List<ReportableVariant> allReportableVariants = Lists.newArrayList();
        for (DriverSomaticVariant somaticDriverVariant : somaticVariantsReport) {
            double adjustedDriverLikelihood = somaticDriverVariant.driverLikelihood();
            for (DriverGermlineVariant germlineVariant : germlineVariantsToReport) {
                if (germlineVariant.variant().gene().equals(somaticDriverVariant.variant().gene())) {
                    adjustedDriverLikelihood = Math.max(adjustedDriverLikelihood, germlineVariant.driverLikelihood());
                }
            }

            allReportableVariants.add(fromSomaticVariant(somaticDriverVariant.variant()).driverLikelihood(adjustedDriverLikelihood)
                    .driverLikelihoodInterpretation(DriverInterpretation.interpret(adjustedDriverLikelihood))
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        for (DriverGermlineVariant driverGermlineVariant : germlineVariantsToReport) {
            double adjustedDriverLikelihood = driverGermlineVariant.driverLikelihood();
            for (DriverSomaticVariant somaticVariant : somaticVariantsReport) {
                if (somaticVariant.variant().gene().equals(driverGermlineVariant.variant().gene())) {
                    adjustedDriverLikelihood = Math.max(adjustedDriverLikelihood, somaticVariant.driverLikelihood());
                }
            }

            allReportableVariants.add(fromGermlineVariant(driverGermlineVariant.variant()).driverLikelihood(adjustedDriverLikelihood)
                    .driverLikelihoodInterpretation(DriverInterpretation.interpret(adjustedDriverLikelihood))
                    .notifyClinicalGeneticist(germlineReportingChoice == LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION
                            && germlineReportingModel.notifyAboutGene(driverGermlineVariant.variant().gene()))
                    .build());
        }

        return allReportableVariants;
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromGermlineVariant(@NotNull ReportableGermlineVariant variant) {
        return ImmutableReportableVariant.builder()
                .gene(variant.gene())
                .position(variant.position())
                .chromosome(variant.chromosome())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalCodingEffect(variant.codingEffect())
                .canonicalHgvsCodingImpact(variant.hgvsCoding())
                .canonicalHgvsProteinImpact(variant.hgvsProtein())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .gDNA(toGDNA(variant.chromosome(), variant.position()))
                .totalCopyNumber(variant.adjustedCopyNumber())
                .alleleCopyNumber(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVaf()))
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromSomaticVariant(@NotNull SomaticVariant variant) {
        return ImmutableReportableVariant.builder()
                .gene(variant.gene())
                .position(variant.position())
                .chromosome(variant.chromosome())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalCodingEffect(variant.canonicalCodingEffect())
                .canonicalHgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .gDNA(toGDNA(variant.chromosome(), variant.position()))
                .totalCopyNumber(variant.adjustedCopyNumber())
                .alleleCopyNumber(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static String toGDNA(@NotNull String chromosome, long position) {
        return chromosome + ":" + position;
    }

    private static double calcAlleleCopyNumber(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
