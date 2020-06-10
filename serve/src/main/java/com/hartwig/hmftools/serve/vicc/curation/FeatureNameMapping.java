package com.hartwig.hmftools.serve.vicc.curation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class FeatureNameMapping {

    @NotNull
    public abstract String originalFeatureName();

    @NotNull
    public abstract String curatedFeatureName();

}
