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
