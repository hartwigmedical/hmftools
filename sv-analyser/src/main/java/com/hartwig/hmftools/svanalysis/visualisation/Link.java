package com.hartwig.hmftools.svanalysis.visualisation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Link {

    public abstract String startChromosome();

    public abstract long startPosition();

    public abstract String endChromosome();

    public abstract long endPosition();
}
