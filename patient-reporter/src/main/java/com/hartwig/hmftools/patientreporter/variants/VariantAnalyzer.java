package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.common.variant.SomaticVariant;

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
    public static VariantAnalyzer of(@NotNull final GeneModel geneModel, @NotNull final MicrosatelliteAnalyzer microsatelliteAnalyzer) {
        return ImmutableVariantAnalyzer.of(new ConsequenceDeterminer(geneModel), microsatelliteAnalyzer);
    }

    @NotNull
    public VariantAnalysis run(@NotNull final List<SomaticVariant> passedVariants) {
        final double indelsPerMb = microsatelliteAnalyzer().analyzeVariants(passedVariants);
        final int mutationalLoad = MutationalLoadAnalyzer.analyzeVariants(passedVariants);

        final ConsequenceOutput consequenceOutput = determiner().run(passedVariants);

        return ImmutableVariantAnalysis.of(passedVariants,
                consequenceOutput.consequentialVariants(),
                consequenceOutput.findings(),
                indelsPerMb,
                mutationalLoad);
    }
}
