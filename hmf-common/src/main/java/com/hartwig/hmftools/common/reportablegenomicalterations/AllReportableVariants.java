package com.hartwig.hmftools.common.reportablegenomicalterations;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bachelor.GermlineReportingModel;
import com.hartwig.hmftools.common.bachelor.GermlineVariant;
import com.hartwig.hmftools.common.drivergene.DriverGeneView;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class AllReportableVariants {

    private AllReportableVariants() {
    }

    @NotNull
    public static List<ReportableVariant> mergeAllSomaticVariants(@NotNull List<SomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull DriverGeneView driverGeneView) {
        List<ReportableVariant> allReportableVariants = Lists.newArrayList();
        for (SomaticVariant variant : somaticVariantsReport) {
            DriverCategory category = driverGeneView.category(variant.gene());
            assert category != null;

            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());
            Double driverLikelihood = null;
            if (catalog != null) {
                driverLikelihood = catalog.driverLikelihood();
            }

            allReportableVariants.add(fromSomaticVariant(variant).driverCategory(category)
                    .driverLikelihood(driverLikelihood)
                    .notifyClinicalGeneticist(false)
                    .build());
        }
        return allReportableVariants;
    }

    @NotNull
    public static List<ReportableVariant> mergeSomaticAndGermlineVariants(@NotNull List<SomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull DriverGeneView driverGeneView,
            @NotNull List<ReportableGermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingChoice germlineReportingChoice) {
        List<ReportableVariant> allReportableVariants = Lists.newArrayList();
        for (SomaticVariant variant : somaticVariantsReport) {
            DriverCategory category = driverGeneView.category(variant.gene());
            assert category != null;

            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());
            Double driverLikelihood = null;
            if (catalog != null) {
                driverLikelihood = catalog.driverLikelihood();
                for (ReportableGermlineVariant germlineVariant : germlineVariantsToReport) {
                    if (germlineVariant.variant().gene().equals(variant.gene())) {
                        driverLikelihood = Math.max(driverLikelihood, germlineVariant.driverLikelihood());
                    }
                }
            }

            allReportableVariants.add(fromSomaticVariant(variant).driverCategory(category)
                    .driverLikelihood(driverLikelihood)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        boolean wantsToBeNotified = germlineReportingChoice == LimsGermlineReportingChoice.ALL
                || germlineReportingChoice == LimsGermlineReportingChoice.ACTIONABLE_ONLY;
        for (ReportableGermlineVariant germlineVariant : germlineVariantsToReport) {
            DriverCategory category = driverGeneView.category(germlineVariant.variant().gene());
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, germlineVariant.variant().gene());
            double driverLikelihood = germlineVariant.driverLikelihood();
            if (catalog != null) {
                driverLikelihood = Math.max(driverLikelihood, catalog.driverLikelihood());
            }
            allReportableVariants.add(fromGermlineVariant(germlineVariant.variant()).driverCategory(category)
                    .driverLikelihood(driverLikelihood)
                    .notifyClinicalGeneticist(wantsToBeNotified && germlineReportingModel.notifyAboutGene(germlineVariant.variant().gene()))
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
                .gDNA(toGDNA(variant))
                .hgvsCodingImpact(variant.hgvsCodingImpact())
                .hgvsProteinImpact(variant.hgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
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
                .gDNA(toGDNA(variant))
                .hgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
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
