package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class ReportableVariantFactory {

    private ReportableVariantFactory() {
    }

    @NotNull
    static List<ReportableVariant> mergeSomaticAndGermlineVariants(@NotNull List<SomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull DriverGeneView driverGeneView,
            @NotNull List<ReportableGermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingChoice) {
        List<ReportableVariant> allReportableVariants = Lists.newArrayList();
        for (SomaticVariant variant : somaticVariantsReport) {
            DriverCategory category = driverGeneView.category(variant.gene());
            assert category != null;

            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());
            Double adjustedDriverLikelihood = null;
            if (catalog != null) {
                adjustedDriverLikelihood = catalog.driverLikelihood();
                for (ReportableGermlineVariant germlineVariant : germlineVariantsToReport) {
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

        for (ReportableGermlineVariant germlineVariant : germlineVariantsToReport) {
            DriverCategory category = driverGeneView.category(germlineVariant.variant().gene());
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

    @NotNull
    private static ImmutableReportableVariant.Builder fromGermlineVariant(@NotNull GermlineVariant variant) {
        return ImmutableReportableVariant.builder()
                .gene(variant.gene())
                .position(variant.position())
                .chromosome(variant.chromosome())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalCodingEffect(variant.codingEffect())
                .canonicalHgvsCodingImpact(variant.hgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.hgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .gDNA(toGDNA(variant))
                .totalPloidy(variant.adjustedCopyNumber())
                .allelePloidy(calcAllelePloidy(variant.adjustedCopyNumber(), variant.adjustedVAF()))
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
                .gDNA(toGDNA(variant))
                .totalPloidy(variant.adjustedCopyNumber())
                .allelePloidy(calcAllelePloidy(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static String toGDNA(@NotNull GenomePosition genomePosition) {
        return genomePosition.chromosome() + ":" + genomePosition.position();
    }

    private static double calcAllelePloidy(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
