package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ReportableVariantAnalyzer {


    @NotNull
    public static List<ReportableVariant> toReportableSomaticVariants(@NotNull List<EnrichedSomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull Map<String, DriverCategory> driverCategoryPerGene,
            @NotNull Set<String> drupActionableGenes, List<GermlineVariant> germlineVariants,
            Map<String,Boolean> germlineGenesReporting) {
        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variants) {
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());

            reportableVariants.add(fromVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryPerGene.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .notifyClinicalGeneticus(false)
                    .build());
        }

        Set<String> notifyGenes = Sets.newHashSet();
        for (String key: germlineGenesReporting.keySet()) {
            if (germlineGenesReporting.get(key)) {
                notifyGenes.add(key);
            }
        }

        for (GermlineVariant germlineVariant : germlineVariants) {
            boolean notify = genesInNotifyClinicalGeneticus(germlineVariants, notifyGenes);

            reportableVariants.add(fromGermline(germlineVariant).isDrupActionable(drupActionableGenes.contains(germlineVariant.gene()))
                    .driverCategory(driverCategoryPerGene.get(germlineVariant.gene()))
                    .driverLikelihood(null)
                    .notifyClinicalGeneticus(notify)
                    .build());

        }
        return reportableVariants;
    }

    private static boolean genesInNotifyClinicalGeneticus(@Nullable List<GermlineVariant> filteredGermlineVariants,
            @NotNull Set<String> notifyGenes) {

        if (filteredGermlineVariants != null) {
            for (GermlineVariant filteredGermlineVariant : filteredGermlineVariants) {
                if (notifyGenes.contains(filteredGermlineVariant.gene())) {
                    return true;
                }
            }
        }
        return false;
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
    private static ImmutableReportableVariant.Builder fromGermline(@NotNull GermlineVariant variantGermline) {
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

    @NotNull
    private static ImmutableReportableVariant.Builder fromVariant(@NotNull EnrichedSomaticVariant variant) {
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

