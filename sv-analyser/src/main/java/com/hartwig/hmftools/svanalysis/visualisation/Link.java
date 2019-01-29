package com.hartwig.hmftools.svanalysis.visualisation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Link {

    enum Type {NORMAL, FOLDBACK, NULL}

    public abstract String sampleId();

    public abstract int clusterId();

    public abstract int chainId();

    public abstract int svId();

    public abstract String startChromosome();

    public abstract long startPosition();

    public abstract int startOrientation();

    public abstract Type startType();

    public abstract String endChromosome();

    public abstract long endPosition();

    public abstract int endOrientation();

    public abstract Type endType();

    public abstract int traverseCount();

}
