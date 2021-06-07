package com.hartwig.hmftools.sage.coverage;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GeneDepth
{

    @NotNull
    String gene();

    double missedVariantLikelihood();

    int[] depthCounts();
}
