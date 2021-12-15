package com.hartwig.hmftools.serve.sources.actin.reader;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActinEntry {

    @NotNull
    public abstract String trial();

    @NotNull
    public abstract ActinRule rule();

    @NotNull
    public abstract List<String> parameters();
}
