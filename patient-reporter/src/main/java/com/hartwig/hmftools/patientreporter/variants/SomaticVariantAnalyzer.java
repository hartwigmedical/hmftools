package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.algo.DriverProbabilityModel;

import org.jetbrains.annotations.NotNull;

public final class SomaticVariantAnalyzer {

    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull final List<EnrichedSomaticVariant> variants, @NotNull Set<String> genePanel,
            @NotNull Map<String, DriverCategory> driverCategoryPerGeneMap) {
        final List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(includeFilter(genePanel, driverCategoryPerGeneMap)).collect(Collectors.toList());
        final double microsatelliteIndelsPerMb = MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants);
        final int tumorMutationalLoad = MutationalLoadAnalyzer.determineTumorMutationalLoad(variants);
        final double tumorMutationalBurden = TumorMutationalBurdenAnalyzer.determineTumorMutationalBurden(variants);

        final List<DriverCatalog> driverCatalog = DriverProbabilityModel.createDriverCatalogForSomaticVariants(variantsToReport);

        return ImmutableSomaticVariantAnalysis.of(variantsToReport,
                driverCatalog,
                microsatelliteIndelsPerMb,
                tumorMutationalLoad,
                tumorMutationalBurden);
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
