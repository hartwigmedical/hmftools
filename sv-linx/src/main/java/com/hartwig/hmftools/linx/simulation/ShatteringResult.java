package com.hartwig.hmftools.linx.simulation;

import org.immutables.value.Value;

@Value.Immutable
public abstract class ShatteringResult
{
    // how many iterations of the specific test were run eg 100K runs simulating 20 breaks
    public abstract int testCount();

    // the run index ie between 0 and the testCount-1
    public abstract int runIndex();

    // the number of segments created by the shattering
    public abstract int segments();

    public abstract int linkedSegments();

    public abstract int exactMatches();

    public abstract int adjacentSegments();
}
