package com.hartwig.hmftools.patientreporter.variants;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.actionability.ActionabilityVariantAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariantAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantAnalyzer.class);

    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull final List<EnrichedSomaticVariant> variants, @NotNull Set<String> genePanel,
            @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap, @NotNull Set<String> drupActionableGenes,
            @Nullable PatientTumorLocation patientTumorLocation, @NotNull ActionabilityAnalyzer actionabilityAnalyzerData)
            throws IOException {
        final double microsatelliteIndelsPerMb = MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants);
        final int tumorMutationalLoad = MutationalLoadAnalyzer.determineTumorMutationalLoad(variants);
        final double tumorMutationalBurden = MutationalBurdenAnalyzer.determineTumorMutationalBurden(variants);

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(DndsDriverGeneLikelihoodSupplier.oncoLikelihood(), variants));
        driverCatalog.addAll(TsgDrivers.drivers(DndsDriverGeneLikelihoodSupplier.tsgLikelihood(), variants));

        final List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(includeFilter(genePanel, driverCategoryPerGeneMap)).collect(Collectors.toList());

        final List<ReportableSomaticVariant> reportableVariants =
                toReportableSomaticVariants(variantsToReport, driverCatalog, driverCategoryPerGeneMap, drupActionableGenes);

        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);
        LOGGER.info("primaryTumorLocation: " + primaryTumorLocation);
        LOGGER.info("cancerTypeMappingReading: " + cancerTypeMappingReading);
        LOGGER.info("doid: " + doidsPrimaryTumorLocation);

        Set<String> actionableGenesVariants = actionabilityAnalyzerData.variantAnalyzer().actionableGenes();

        LOGGER.info("Looking for actionability for somatic variants");
        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant = ActionabilityVariantAnalyzer.findEvidenceForVariants(
                actionableGenesVariants,
                variants,
                doidsPrimaryTumorLocation,
                actionabilityAnalyzerData);


        return ImmutableSomaticVariantAnalysis.of(reportableVariants,
                evidencePerVariant,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden);
    }

    @NotNull
    private static List<ReportableSomaticVariant> toReportableSomaticVariants(@NotNull List<EnrichedSomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull Map<String, DriverCategory> driverCategoryPerGene,
            @NotNull Set<String> drupActionableGenes) {
        List<ReportableSomaticVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variants) {
            DriverCatalog catalog = catalogEntryForVariant(driverCatalog, variant);

            reportableVariants.add(fromVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryPerGene.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .build());
        }

        return reportableVariants;
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
                // KODU: If a variant has uncertain driver category we should always report.
                return TSG_CODING_EFFECTS_TO_REPORT.contains(effect) || ONCO_CODING_EFFECTS_TO_REPORT.contains(effect);
            }
        };
    }
}
