package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.HmfReporterData;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
public abstract class SomaticVariantAnalyzer {

    @NotNull
    abstract ConsequenceDeterminer determiner();

    @NotNull
    abstract MicrosatelliteAnalyzer microsatelliteAnalyzer();

    @NotNull
    public static SomaticVariantAnalyzer of(@NotNull HmfReporterData reporterData) {
        final Set<String> transcriptsToInclude = reporterData.panelGeneModel().transcriptMap().keySet();
        return of(transcriptsToInclude, ImmutableMicrosatelliteAnalyzer.of(reporterData.refGenomeFastaFile()));
    }

    @VisibleForTesting
    @NotNull
    public static SomaticVariantAnalyzer of(@NotNull Set<String> transcripts, @NotNull MicrosatelliteAnalyzer microsatelliteAnalyzer) {
        return ImmutableSomaticVariantAnalyzer.of(new ConsequenceDeterminer(transcripts), microsatelliteAnalyzer);
    }

    @NotNull
    public SomaticVariantAnalysis run(@NotNull final List<EnrichedSomaticVariant> variants) {
        final double indelsPerMb = microsatelliteAnalyzer().analyzeVariants(variants);
        final int mutationalLoad = MutationalLoadAnalyzer.determineMutationalLoad(variants);

        final List<VariantReport> variantReports = determiner().run(variants);

        return ImmutableSomaticVariantAnalysis.of(variantReports, indelsPerMb, mutationalLoad);
    }
}
