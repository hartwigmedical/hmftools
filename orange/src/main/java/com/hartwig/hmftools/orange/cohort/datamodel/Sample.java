package com.hartwig.hmftools.orange.cohort.datamodel;

import java.util.Set;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Sample {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract Set<String> doids();

}
