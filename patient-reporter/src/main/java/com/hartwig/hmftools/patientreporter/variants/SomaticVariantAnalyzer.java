package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class SomaticVariantAnalyzer {

    private static final List<CodingEffect> CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull final List<EnrichedSomaticVariant> variants, @NotNull Set<String> genePanel) {
        final List<EnrichedSomaticVariant> variantsToReport =
                variants.stream().filter(variantFilter(genePanel)).collect(Collectors.toList());
        final double indelsPerMb = MicrosatelliteAnalyzer.determineMicrosatelliteIndelsPerMb(variants);
        final int mutationalLoad = MutationalLoadAnalyzer.determineMutationalLoad(variants);

        return ImmutableSomaticVariantAnalysis.of(variantsToReport, indelsPerMb, mutationalLoad);
    }

    @NotNull
    private static Predicate<EnrichedSomaticVariant> variantFilter(@NotNull Set<String> genePanel) {
        return variant -> CODING_EFFECTS_TO_REPORT.contains(variant.canonicalCodingEffect()) && genePanel.contains(variant.gene())
                && !variant.isFiltered();
    }
}
