package com.hartwig.hmftools.orange.cohort.mapping;

import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CohortMapping {

    @NotNull
    public abstract String cancerType();

    public abstract int preferenceRank();

    @NotNull
    public abstract MappingRule rule();

    @NotNull
    public abstract Set<String> include();

    @NotNull
    public abstract Set<String> exclude();
}
