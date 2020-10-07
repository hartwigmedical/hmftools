package com.hartwig.hmftools.serve.vicc.annotation;

import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneRangeAnnotation {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String chromosome();

    public abstract long start();

    public abstract long end();

    @NotNull
    public abstract MutationTypeFilter mutationType();

}
