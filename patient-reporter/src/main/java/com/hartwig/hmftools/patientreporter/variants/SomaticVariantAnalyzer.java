package com.hartwig.hmftools.patientreporter.variants;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeMappingReading;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVs;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsEvidenceItems;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRange;
import com.hartwig.hmftools.common.actionability.somaticvariant.ActionabilityRangeEvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.EvidenceItem;
import com.hartwig.hmftools.common.actionability.somaticvariant.VariantEvidenceItems;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;

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
            @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap, @Nullable PatientTumorLocation patientTumorLocation,
            @NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull List<GeneFusion> fusions,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) throws IOException {
        final List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(includeFilter(genePanel, driverCategoryPerGeneMap)).collect(Collectors.toList());
        final double microsatelliteIndelsPerMb = MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants);
        final int tumorMutationalLoad = MutationalLoadAnalyzer.determineTumorMutationalLoad(variants);
        final double tumorMutationalBurden = MutationalBurdenAnalyzer.determineTumorMutationalBurden(variants);

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(DndsDriverGeneLikelihoodSupplier.oncoLikelihood(), variants));
        driverCatalog.addAll(TsgDrivers.drivers(DndsDriverGeneLikelihoodSupplier.tsgLikelihood(), variants));

        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
        CancerTypeMappingReading cancerTypeMappingReading = CancerTypeMappingReading.readingFile();
        String doidsPrimaryTumorLocation = cancerTypeMappingReading.doidsForPrimaryTumorLocation(primaryTumorLocation);
        LOGGER.info("primaryTumorLocation: " + primaryTumorLocation);
        LOGGER.info("cancerTypeMappingReading: " + cancerTypeMappingReading);
        LOGGER.info("doid: " + doidsPrimaryTumorLocation);

        Set<String> actionableGenesVariants = actionabilityAnalyzerData.variantAnalyzer().actionableGenes();
        Set<String> actionableGenesCNVS = actionabilityAnalyzerData.cnvAnalyzer().actionableGenes();

        LOGGER.info("evidence items variants");
        final List<EvidenceItem> variant = Lists.newArrayList();
        Map<EnrichedSomaticVariant, VariantEvidenceItems> evidencePerVariant = ActionabilityVariantAnalyzer.detectVariants(
                actionableGenesVariants,
                variants,
                doidsPrimaryTumorLocation,
                actionabilityAnalyzerData);

        for (Map.Entry<EnrichedSomaticVariant, VariantEvidenceItems> entry : evidencePerVariant.entrySet()) {
            variant.addAll(entry.getValue().onLabel());
            variant.addAll(entry.getValue().offLabel());
        }

        LOGGER.info("evidence items variants ranges");
        final List<ActionabilityRange> variantRange = Lists.newArrayList();
        Map<EnrichedSomaticVariant, ActionabilityRangeEvidenceItem> evidencePerVariantRanges =
                ActionabilityVariantAnalyzer.detectVariantsRanges(actionableGenesVariants,
                        variants,
                        doidsPrimaryTumorLocation,
                        actionabilityAnalyzerData);

        for (Map.Entry<EnrichedSomaticVariant, ActionabilityRangeEvidenceItem> entryRange : evidencePerVariantRanges.entrySet()) {
            variantRange.addAll(entryRange.getValue().onLabel());
            variantRange.addAll(entryRange.getValue().offLabel());
        }

        LOGGER.info("evidence items CNVs");
        final List<ActionabilityCNVs> CNVs = Lists.newArrayList();
        Map<GeneCopyNumber, ActionabilityCNVsEvidenceItems> evidencePerVariantCNVs = ActionabilityVariantAnalyzer.detectCNVs(
                actionableGenesCNVS,
                geneCopyNumbers,
                doidsPrimaryTumorLocation,
                actionabilityAnalyzerData);

        for (Map.Entry<GeneCopyNumber, ActionabilityCNVsEvidenceItems> entryCNVs : evidencePerVariantCNVs.entrySet()) {
            CNVs.addAll(entryCNVs.getValue().onLabel());
            CNVs.addAll(entryCNVs.getValue().offLabel());
        }

        return ImmutableSomaticVariantAnalysis.of(variantsToReport,
                driverCatalog,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden,
                variant,
                variantRange,
                CNVs,
                evidencePerVariant,
                evidencePerVariantRanges,
                evidencePerVariantCNVs);
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
