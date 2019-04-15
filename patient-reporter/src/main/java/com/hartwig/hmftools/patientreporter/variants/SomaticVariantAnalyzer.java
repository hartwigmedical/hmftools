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
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.germline.ImmutableGermlineVariant;

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
            List<GermlineVariant> filteredGermlineVariant, @NotNull Set<String> notifyGeneticus) {
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

        List<ReportableSomaticVariant> reportableVariants = toReportableSomaticVariants(variantsToReport,
                driverCatalog,
                driverCategoryPerGeneMap,
                drupActionableGenes,
                filteredGermlineVariant, notifyGeneticus);

        return ImmutableSomaticVariantAnalysis.of(reportableVariants,
                filteredEvidence,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden);
    }

    @NotNull
    private static List<ReportableSomaticVariant> toReportableSomaticVariants(@NotNull List<EnrichedSomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull Map<String, DriverCategory> driverCategoryPerGene,
            @NotNull Set<String> drupActionableGenes, List<GermlineVariant> filteredGermlineVariant, @NotNull Set<String> notifyGeneticus) {
        List<ReportableSomaticVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variants) {
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant);
            List<GermlineVariant> germlineVariant = fromGermlineReporting(filteredGermlineVariant, variant);
            boolean notify = genesInNotifyClinicalGeneticus(filteredGermlineVariant, notifyGeneticus);

            reportableVariants.add(fromVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryPerGene.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .germlineVariant(germlineVariant)
                    .SomaticOrGermline(germlineVariant.size() == 0 ? "somatic" : "germline")
                    .notifyClinicalGeneticus(notify)
                    .build());
        }

        return reportableVariants;
    }

    private static boolean genesInNotifyClinicalGeneticus(@Nullable List<GermlineVariant> filteredGermlineVariants, @NotNull Set<String> notifyGeneticus) {
        if (filteredGermlineVariants != null) {
            for (GermlineVariant filteredGermlineVariant : filteredGermlineVariants) {
                if (notifyGeneticus.contains(filteredGermlineVariant.gene())) {
                    return true;
                }
            }
        }
        return false;
    }

    @NotNull
    private static List<GermlineVariant> fromGermlineReporting(@Nullable List<GermlineVariant> filteredGermlineVariants,
            @NotNull SomaticVariant variant) {
        List<GermlineVariant> germlineVariant = Lists.newArrayList();
        if (filteredGermlineVariants != null) {
            for (GermlineVariant filteredGermlineVariant : filteredGermlineVariants) {

                return Lists.newArrayList(ImmutableGermlineVariant.builder()
                        .passFilter(filteredGermlineVariant.passFilter())
                        .gene(filteredGermlineVariant.gene())
                        .hgvsCodingImpact(filteredGermlineVariant.hgvsCodingImpact())
                        .hgvsProteinImpact(filteredGermlineVariant.hgvsProteinImpact())
                        .totalReadCount(filteredGermlineVariant.totalReadCount())
                        .alleleReadCount(filteredGermlineVariant.alleleReadCount())
                        .germlineStatus(filteredGermlineVariant.germlineStatus())
                        .adjustedVAF(filteredGermlineVariant.adjustedVAF())
                        .adjustedCopyNumber(filteredGermlineVariant.adjustedCopyNumber())
                        .minorAllelePloidy(filteredGermlineVariant.minorAllelePloidy())
                        .biallelic(filteredGermlineVariant.biallelic())
                        .build());
            }
        }
        return germlineVariant;
    }

    @Nullable
    private static DriverCatalog catalogEntryForVariant(@NotNull List<DriverCatalog> driverCatalogList, @NotNull SomaticVariant variant) {
        for (DriverCatalog entry : driverCatalogList) {
            if (entry.gene().equals(variant.gene())) {
                return entry;
            }
        }
        return null;
    }

    @NotNull
    private static ImmutableReportableSomaticVariant.Builder fromVariant(@NotNull EnrichedSomaticVariant variant) {
        return ImmutableReportableSomaticVariant.builder()
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
