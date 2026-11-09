package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// A probe segment consisting of a reference genome region. If orientation is REVERSE the bases are reverse complemented.
public record RefSegment(ChrBaseRegion region, Orientation orientation) implements SequenceSegment
{
    public RefSegment
    {
        requireNonNull(region);
        requireNonNull(orientation);
        if(!region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region: " + region);
        }
    }

    @Override
    public int baseLength()
    {
        return region.baseLength();
    }

    // Whether this segment is immediately followed in probe-sequence space by the given segment along the genome with the same orientation,
    // i.e. the two are directly adjacent and could be merged into a single region.
    public boolean isGenomeAdjacentTo(final RefSegment next)
    {
        if(orientation != next.orientation || !region.chromosome().equals(next.region.chromosome()))
        {
            return false;
        }
        // For FORWARD, the sequence runs low->high, so adjacency means this region ends just before the next begins.
        // For REVERSE, the sequence runs high->low, so adjacency means the next region ends just before this one begins.
        return orientation == Orientation.FORWARD
                ? region.end() + 1 == next.region.start()
                : next.region.end() + 1 == region.start();
    }

    @Override
    public int typeRank()
    {
        return 0;
    }

    @Override
    public int compareSameType(final SequenceSegment other)
    {
        RefSegment o = (RefSegment) other;
        int compare = region.compareTo(o.region);
        return compare != 0 ? compare : orientation.compareTo(o.orientation);
    }
}
