package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class ReportableVariantFactory {

    private ReportableVariantFactory() {
    }

    @NotNull
    static List<ReportableVariant> mergeSomaticAndGermlineVariants(@NotNull List<SomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull List<DriverGene> driverGenes,
            @NotNull List<ReportableGermlineVariantExtended> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingChoice) {
        List<ReportableVariant> allReportableVariants = Lists.newArrayList();
        for (SomaticVariant variant : somaticVariantsReport) {
            DriverCategory category = catalogDriverGene(driverGenes, variant.gene());

            assert category != null;

            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());
            Double adjustedDriverLikelihood = null;
            if (catalog != null) {
                adjustedDriverLikelihood = catalog.driverLikelihood();
                for (ReportableGermlineVariantExtended germlineVariant : germlineVariantsToReport) {
                    if (germlineVariant.variant().gene().equals(variant.gene())) {
                        adjustedDriverLikelihood = Math.max(adjustedDriverLikelihood, germlineVariant.driverLikelihood());
                    }
                }
            }

            allReportableVariants.add(fromSomaticVariant(variant).driverCategory(category)
                    .driverLikelihood(adjustedDriverLikelihood)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        for (ReportableGermlineVariantExtended germlineVariant : germlineVariantsToReport) {
            DriverCategory category = catalogDriverGene(driverGenes, germlineVariant.variant().gene());
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, germlineVariant.variant().gene());
            double adjustedDriverLikelihood = germlineVariant.driverLikelihood();
            if (catalog != null) {
                adjustedDriverLikelihood = Math.max(adjustedDriverLikelihood, catalog.driverLikelihood());
            }
            allReportableVariants.add(fromGermlineVariant(germlineVariant.variant()).driverCategory(category)
                    .driverLikelihood(adjustedDriverLikelihood)
                    .notifyClinicalGeneticist(germlineReportingChoice == LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION
                            && germlineReportingModel.notifyAboutGene(germlineVariant.variant().gene()))
                    .build());
        }

        return allReportableVariants;
    }

    @Nullable
    private static DriverCatalog catalogEntryForVariant(@NotNull List<DriverCatalog> driverCatalogList, @NotNull String gene) {
        for (DriverCatalog entry : driverCatalogList) {
            if (entry.gene().equals(gene)) {
                return entry;
            }
        }
        return null;
    }

    @Nullable
    private static DriverCategory catalogDriverGene(@NotNull List<DriverGene> driverGenes, @NotNull String gene) {
        for (DriverGene driverGene: driverGenes) {
            if (driverGene.gene().equals(gene)) {
                return driverGene.likelihoodType();
            }
        }
        return null;
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
                .totalPloidy(variant.adjustedCopyNumber())
                .allelePloidy(calcAllelePloidy(variant.adjustedCopyNumber(), variant.adjustedVaf()))
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
                .totalPloidy(variant.adjustedCopyNumber())
                .allelePloidy(calcAllelePloidy(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static String toGDNA(@NotNull String chromosome, long position) {
        return chromosome + ":" + position;
    }

    private static double calcAllelePloidy(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
