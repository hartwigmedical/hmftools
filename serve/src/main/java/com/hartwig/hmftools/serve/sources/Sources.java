package com.hartwig.hmftools.serve.sources;

import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class Sources {

    @NotNull
    public abstract String sourceEvent();

    @NotNull
    public abstract Knowledgebase source();
}
