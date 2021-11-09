package com.hartwig.hmftools.orange.cohort;

import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CohortMapping {

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract CohortMappingRule mappingRule();

    @NotNull
    public abstract Set<String> include();

    @NotNull
    public abstract Set<String> exclude();
}
