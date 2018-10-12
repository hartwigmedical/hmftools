package com.hartwig.hmftools.patientreporter.variants;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVs;
import com.hartwig.hmftools.common.actionability.cnv.ActionabilityCNVsEvidenceItems;
import com.hartwig.hmftools.common.actionability.fusion.ActionabilityFusionPairs;
import com.hartwig.hmftools.common.actionability.fusion.FusionEvidenceItems;
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
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.copynumber.PurpleAnalysis;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariantAnalyzer {
    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull final List<EnrichedSomaticVariant> variants, @NotNull Set<String> genePanel,
            @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap, @Nullable PatientTumorLocation patientTumorLocation,
            @NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull List<GeneFusion> fusions) throws IOException {
        final List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(includeFilter(genePanel, driverCategoryPerGeneMap)).collect(Collectors.toList());
        final double microsatelliteIndelsPerMb = MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants);
        final int tumorMutationalLoad = MutationalLoadAnalyzer.determineTumorMutationalLoad(variants);
        final double tumorMutationalBurden = TumorMutationalBurdenAnalyzer.determineTumorMutationalBurden(variants);

        final List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(DndsDriverGeneLikelihoodSupplier.oncoLikelihood(), variants));
        driverCatalog.addAll(TsgDrivers.drivers(DndsDriverGeneLikelihoodSupplier.tsgLikelihood(), variants));

        final List<EvidenceItem> variant = Lists.newArrayList();
        Map<EnrichedSomaticVariant, VariantEvidenceItems> evidencePerVariant =
                ActionabilityVariantAnalyzer.detectVariants(variants, patientTumorLocation);

        final List<ActionabilityRange> variantRange = Lists.newArrayList();
        Map<EnrichedSomaticVariant, ActionabilityRangeEvidenceItem> evidencePerVariantRanges =
                ActionabilityVariantAnalyzer.detectVariantsRanges(variants, patientTumorLocation);

        final List<ActionabilityCNVs> CNVs = Lists.newArrayList();
        Map<GeneCopyNumber, ActionabilityCNVsEvidenceItems> evidencePerVariantCNVs =
                ActionabilityVariantAnalyzer.detectCNVs(geneCopyNumbers, patientTumorLocation);

        for (Map.Entry<EnrichedSomaticVariant, VariantEvidenceItems> entry : evidencePerVariant.entrySet()) {
            variant.addAll(entry.getValue().onLabel());
            variant.addAll(entry.getValue().offLabel());
        }

        for (Map.Entry<EnrichedSomaticVariant, ActionabilityRangeEvidenceItem> entryRange : evidencePerVariantRanges.entrySet()) {
            variantRange.addAll(entryRange.getValue().onLabel());
            variantRange.addAll(entryRange.getValue().offLabel());
        }

        for (Map.Entry<GeneCopyNumber, ActionabilityCNVsEvidenceItems> entryRange : evidencePerVariantCNVs.entrySet()) {
            CNVs.addAll(entryRange.getValue().onLabel());
            CNVs.addAll(entryRange.getValue().offLabel());
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
