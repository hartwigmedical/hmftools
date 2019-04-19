package com.hartwig.hmftools.patientreporter.variants;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.germline.FilterGermlineVariants;
import com.hartwig.hmftools.patientreporter.germline.GermlineGenesReporting;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariantAnalyzer {
    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);
    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantAnalyzer.class);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull List<EnrichedSomaticVariant> variants, @NotNull Set<String> genePanel,
            @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap, @NotNull Set<String> drupActionableGenes,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzer, @Nullable PatientTumorLocation patientTumorLocation,
            List<GermlineVariant> germlineVariants, @NotNull Set<String> notifyGeneticus, GermlineGenesReporting germlineGenesReporting,
            @NotNull String sample, @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        double microsatelliteIndelsPerMb = MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants);
        int tumorMutationalLoad = MutationalLoadAnalyzer.determineTumorMutationalLoad(variants);
        double tumorMutationalBurden = MutationalBurdenAnalyzer.determineTumorMutationalBurden(variants);

        List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(DndsDriverGeneLikelihoodSupplier.oncoLikelihood(), variants));
        driverCatalog.addAll(TsgDrivers.drivers(DndsDriverGeneLikelihoodSupplier.tsgLikelihood(), variants));

        List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(includeFilter(genePanel, driverCategoryPerGeneMap)).collect(Collectors.toList());

        String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : null;
        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant =
                actionabilityAnalyzer.evidenceForSomaticVariants(variants, primaryTumorLocation);

        List<EvidenceItem> filteredEvidence = ReportableEvidenceItemFactory.reportableFlatList(evidencePerVariant);

        // Add all variants with filtered evidence that have not previously been added.
        for (Map.Entry<EnrichedSomaticVariant, List<EvidenceItem>> entry : evidencePerVariant.entrySet()) {
            EnrichedSomaticVariant variant = entry.getKey();
            if (!variantsToReport.contains(variant) && !Collections.disjoint(entry.getValue(), filteredEvidence)) {
                variantsToReport.add(variant);
            }
        }

        List<ReportableVariant> reportableVariants = toReportableSomaticVariants(variantsToReport,
                driverCatalog,
                driverCategoryPerGeneMap,
                drupActionableGenes,
                germlineVariants,
                notifyGeneticus,
                germlineGenesReporting,
                geneCopyNumbers,
                sample);

        return ImmutableSomaticVariantAnalysis.of(reportableVariants,
                filteredEvidence,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden);
    }

    @NotNull
    private static List<ReportableVariant> toReportableSomaticVariants(@NotNull List<EnrichedSomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull Map<String, DriverCategory> driverCategoryPerGene,
            @NotNull Set<String> drupActionableGenes, List<GermlineVariant> germlineVariants, @NotNull Set<String> notifyGeneticus,
            GermlineGenesReporting germlineGenesReporting,
            @NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull String sampleId) {
        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variants) {
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant.gene());
            boolean notify = genesInNotifyClinicalGeneticus(germlineVariants, notifyGeneticus);

            reportableVariants.add(fromVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryPerGene.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .somaticOrGermline("somatic")
                    .notifyClinicalGeneticus(notify)
                    .build());
        }

        final List<GermlineVariant> filteredGermlineVariants = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariants,
                germlineGenesReporting,
                driverCategoryPerGene,
                geneCopyNumbers,
                sampleId,
                variants);

        for (GermlineVariant germlineVariant : filteredGermlineVariants) {
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, germlineVariant.gene());
            boolean notify = genesInNotifyClinicalGeneticus(filteredGermlineVariants, notifyGeneticus);

            reportableVariants.add(fromGermline(germlineVariant).isDrupActionable(drupActionableGenes.contains(germlineVariant.gene()))
                    .driverCategory(driverCategoryPerGene.get(germlineVariant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .somaticOrGermline("germline")
                    .notifyClinicalGeneticus(notify)
                    .build());

        }
        return reportableVariants;
    }

    private static boolean genesInNotifyClinicalGeneticus(@Nullable List<GermlineVariant> filteredGermlineVariants,
            @NotNull Set<String> notifyGeneticus) {
        if (filteredGermlineVariants != null) {
            for (GermlineVariant filteredGermlineVariant : filteredGermlineVariants) {
                if (notifyGeneticus.contains(filteredGermlineVariant.gene())) {
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

    @NotNull
    private static Predicate<EnrichedSomaticVariant> includeFilter(@NotNull Set<String> genePanel,
            @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap) {
        return variant -> {
            if (variant.isFiltered()) {
                return false;
            }

            if (!genePanel.contains(variant.gene())) {
                return false;
            }

            CodingEffect effect = variant.canonicalCodingEffect();
            if (driverCategoryPerGeneMap.get(variant.gene()) == DriverCategory.TSG) {
                return TSG_CODING_EFFECTS_TO_REPORT.contains(effect);
            } else if (driverCategoryPerGeneMap.get(variant.gene()) == DriverCategory.ONCO) {
                return ONCO_CODING_EFFECTS_TO_REPORT.contains(effect);
            } else {
                // If a variant has uncertain driver category we should always report.
                return TSG_CODING_EFFECTS_TO_REPORT.contains(effect) || ONCO_CODING_EFFECTS_TO_REPORT.contains(effect);
            }
        };
    }
}
