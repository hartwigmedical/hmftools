package com.hartwig.hmftools.svtools.simulation;

import org.immutables.value.Value;

@Value.Immutable
public abstract class ShatteringResult
{
    // the run index ie between 0 and the testCount-1
    public abstract int runIndex();

    // the number of segments created by the shattering
    public abstract int segments();

    public abstract int linkedSegments();

    // number of breaks repaired exactly back to original sequence
    public abstract int exactRepairs();

    // number of segments in the right order even if flipped around
    public abstract int adjacentSegments();

    // after factoring out exact repairs, number of distinct segments (ie chain links)
    public abstract int inferredLinks();

    // segments lost after counting contiguous lost sections as a single unit
    public abstract int inferredLost();

    // string representation of the links made between segments
    public abstract String linkStr();

    public boolean equals(final ShatteringResult other)
    {
        return segments() == other.segments() && linkedSegments() == other.linkedSegments()
                && exactRepairs() == other.exactRepairs() && adjacentSegments() == other.adjacentSegments()
                && inferredLinks() == other.inferredLinks() && inferredLost() == other.inferredLost();
    }

}
