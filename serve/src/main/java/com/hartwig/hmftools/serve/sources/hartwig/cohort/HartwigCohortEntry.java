package com.hartwig.hmftools.serve.sources.hartwig.cohort;

import com.hartwig.hmftools.serve.sources.hartwig.HartwigEntry;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HartwigCohortEntry implements HartwigEntry {

    @NotNull
    public abstract String chromosome();

    public abstract long position();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String proteinAnnotation();
}
