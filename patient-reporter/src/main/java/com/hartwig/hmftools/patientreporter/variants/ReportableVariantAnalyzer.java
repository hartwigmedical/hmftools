package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableVariantAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(ReportableVariantAnalyzer.class);

    private ReportableVariantAnalyzer() {
    }

    @NotNull
    public static List<ReportableVariant> mergeSomaticAndGermlineVariants(@NotNull List<EnrichedSomaticVariant> somaticVariantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull DriverGeneView driverGeneView, @NotNull Set<String> drupActionableGenes,
            @NotNull List<GermlineVariant> germlineVariantsToReport, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingChoice germlineReportingChoice) {
        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : somaticVariantsReport) {
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());
            if (catalog == null) {
                LOGGER.warn("No driver entry found for gene {}!", variant.gene());
            }

            reportableVariants.add(fromSomaticVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverGeneView.category(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        boolean wantsToBeNotified = germlineReportingChoice == LimsGermlineReportingChoice.ALL
                || germlineReportingChoice == LimsGermlineReportingChoice.ACTIONABLE_ONLY;
        for (GermlineVariant germlineVariant : germlineVariantsToReport) {
            reportableVariants.add(fromGermlineVariant(germlineVariant).isDrupActionable(drupActionableGenes.contains(germlineVariant.gene()))
                    .driverCategory(driverGeneView.category(germlineVariant.gene()))
                    .driverLikelihood(null)
                    .notifyClinicalGeneticist(wantsToBeNotified && germlineReportingModel.notifyAboutGene(germlineVariant.gene()))
                    .build());

        }
        return reportableVariants;
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
                .gDNA(toGDNA(variant))
                .hgvsCodingImpact(variant.hgvsCodingImpact())
                .hgvsProteinImpact(variant.hgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonalLikelihood(1D)
                .adjustedCopyNumber(variant.adjustedCopyNumber())
                .adjustedVAF(variant.adjustedVAF())
                .minorAllelePloidy(variant.minorAllelePloidy())
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromSomaticVariant(@NotNull EnrichedSomaticVariant variant) {
        return ImmutableReportableVariant.builder()
                .gene(variant.gene())
                .gDNA(toGDNA(variant))
                .hgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .adjustedCopyNumber(variant.adjustedCopyNumber())
                .adjustedVAF(variant.adjustedVAF())
                .minorAllelePloidy(variant.minorAllelePloidy())
                .biallelic(variant.biallelic());
    }

    @NotNull
    private static String toGDNA(@NotNull GenomePosition genomePosition) {
        return genomePosition.chromosome() + ":" + genomePosition.position();
    }
}

