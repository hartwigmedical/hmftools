package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.HmfReporterData;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
public abstract class VariantAnalyzer {

    @NotNull
    protected abstract ConsequenceDeterminer determiner();

    @NotNull
    protected abstract MicrosatelliteAnalyzer microsatelliteAnalyzer();

    @NotNull
    public static VariantAnalyzer of(@NotNull HmfReporterData reporterData) {
        Set<String> transcriptsToInclude = reporterData.panelGeneModel().transcriptMap().keySet();
        return of(transcriptsToInclude, reporterData.microsatelliteAnalyzer());
    }

    @VisibleForTesting
    @NotNull
    public static VariantAnalyzer of(@NotNull Set<String> transcripts, @NotNull MicrosatelliteAnalyzer microsatelliteAnalyzer) {
        return ImmutableVariantAnalyzer.of(new ConsequenceDeterminer(transcripts), microsatelliteAnalyzer);
    }

    @NotNull
    public VariantAnalysis run(@NotNull final List<SomaticVariant> passedVariants) {
        final double indelsPerMb = microsatelliteAnalyzer().analyzeVariants(passedVariants);
        final int mutationalLoad = MutationalLoadAnalyzer.analyzeVariants(passedVariants);

        final List<VariantReport> variantReports = determiner().run(passedVariants);

        return ImmutableVariantAnalysis.of(passedVariants, variantReports, indelsPerMb, mutationalLoad);
    }
}
