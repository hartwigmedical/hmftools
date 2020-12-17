package com.hartwig.hmftools.serve.sources.vicc.curation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FeatureCurationValues {

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract String featureName();

}
