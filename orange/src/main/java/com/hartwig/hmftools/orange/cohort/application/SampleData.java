package com.hartwig.hmftools.orange.cohort.application;

import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class SampleData {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract Set<String> doids();

}
