package com.hartwig.hmftools.svannotation.annotations;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneDisruption {

    public abstract boolean reportable();

    public abstract Transcript linkedAnnotation();
}
