package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableVariantAnalyzer {

    private ReportableVariantAnalyzer() {
    }

    @NotNull
    public static List<ReportableVariant> mergeSomaticAndGermlineVariants(@NotNull List<EnrichedSomaticVariant> variantsReport,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull Map<String, DriverCategory> driverCategoryPerGene,
            @NotNull Set<String> drupActionableGenes, List<GermlineVariant> filteredGermlineVariants,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull LimsGermlineReportingChoice germlineReportingChoice) {
        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variantsReport) {
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());

            reportableVariants.add(fromSomaticVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryPerGene.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        boolean wantsToBeNotified = germlineReportingChoice == LimsGermlineReportingChoice.ALL
                || germlineReportingChoice == LimsGermlineReportingChoice.ACTIONABLE_ONLY;
        for (GermlineVariant germlineVariant : filteredGermlineVariants) {
            reportableVariants.add(fromGermlineVariant(germlineVariant).isDrupActionable(drupActionableGenes.contains(germlineVariant.gene()))
                    .driverCategory(driverCategoryPerGene.get(germlineVariant.gene()))
                    .driverLikelihood(null)
                    .notifyClinicalGeneticist(
                            wantsToBeNotified && germlineReportingModel.notifyAboutGene(germlineVariant.gene()))
                    .build());

        }
        return reportableVariants;
    }

    @VisibleForTesting
    @Nullable
    public static DriverCatalog catalogEntryForVariant(@NotNull List<DriverCatalog> driverCatalogList, @NotNull String gene) {
        for (DriverCatalog entry : driverCatalogList) {
            if (entry.gene().equals(gene)) {
                return entry;
            }
        }
        return null;
    }

    @VisibleForTesting
    @NotNull
    public static ImmutableReportableVariant.Builder fromGermlineVariant(@NotNull GermlineVariant variantGermline) {
        return ImmutableReportableVariant.builder()
                .gene(variantGermline.gene())
                .hgvsCodingImpact(variantGermline.hgvsCodingImpact())
                .hgvsProteinImpact(variantGermline.hgvsProteinImpact())
                .totalReadCount(variantGermline.totalReadCount())
                .alleleReadCount(variantGermline.alleleReadCount())
                .hotspot(Hotspot.NON_HOTSPOT)
                .clonality(Clonality.CLONAL)
                .adjustedCopyNumber(variantGermline.adjustedCopyNumber())
                .adjustedVAF(variantGermline.adjustedVAF())
                .minorAllelePloidy(variantGermline.minorAllelePloidy())
                .biallelic(variantGermline.biallelic());
    }

    @VisibleForTesting
    @NotNull
    public static ImmutableReportableVariant.Builder fromSomaticVariant(@NotNull EnrichedSomaticVariant variant) {
        return ImmutableReportableVariant.builder()
                .gene(variant.gene())
                .hgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .hgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .hotspot(variant.hotspot())
                .clonality(variant.clonality())
                .adjustedCopyNumber(variant.adjustedCopyNumber())
                .adjustedVAF(variant.adjustedVAF())
                .minorAllelePloidy(variant.minorAllelePloidy())
                .biallelic(variant.biallelic());
    }
}

