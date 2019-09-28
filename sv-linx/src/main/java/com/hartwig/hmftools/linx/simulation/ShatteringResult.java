package com.hartwig.hmftools.linx.simulation;

import org.immutables.value.Value;

@Value.Immutable
public abstract class ShatteringResult
{
    public abstract int testCount();

    public abstract int runIndex();

    public abstract int segments();

    public abstract int linkedSegments();

    public abstract int exactMatches();

    public abstract int adjacentSegments();
}
