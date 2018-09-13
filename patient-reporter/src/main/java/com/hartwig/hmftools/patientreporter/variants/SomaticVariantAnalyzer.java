package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.HmfReporterData;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
public abstract class SomaticVariantAnalyzer {

    private static final List<CodingEffect> CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    @NotNull
    abstract Set<String> genePanel();

    @NotNull
    abstract MicrosatelliteAnalyzer microsatelliteAnalyzer();

    @NotNull
    public static SomaticVariantAnalyzer of(@NotNull HmfReporterData reporterData) {
        final Set<String> transcriptsToInclude = reporterData.panelGeneModel().panel();
        return of(transcriptsToInclude, ImmutableMicrosatelliteAnalyzer.of(reporterData.refGenomeFastaFile()));
    }

    @VisibleForTesting
    @NotNull
    static SomaticVariantAnalyzer of(@NotNull Set<String> genePanel, @NotNull MicrosatelliteAnalyzer microsatelliteAnalyzer) {
        return ImmutableSomaticVariantAnalyzer.of(genePanel, microsatelliteAnalyzer);
    }

    @NotNull
    public SomaticVariantAnalysis run(@NotNull final List<EnrichedSomaticVariant> variants) {
        final List<EnrichedSomaticVariant> variantsToReport = variants.stream().filter(variantFilter()).collect(Collectors.toList());
        final double indelsPerMb = microsatelliteAnalyzer().analyzeVariants(variants);
        final int mutationalLoad = MutationalLoadAnalyzer.determineMutationalLoad(variants);

        return ImmutableSomaticVariantAnalysis.of(variantsToReport, indelsPerMb, mutationalLoad);
    }

    @NotNull
    private Predicate<EnrichedSomaticVariant> variantFilter() {
        return variant -> CODING_EFFECTS_TO_REPORT.contains(variant.canonicalCodingEffect()) && genePanel().contains(variant.gene())
                && !variant.isFiltered();
    }
}
