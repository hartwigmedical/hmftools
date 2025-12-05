package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionIntersection;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Range within a probe sequence which the probe is intending to target.
public record TargetedRange(
        // 0-indexed to start of the probe sequence. Inclusive.
        int startOffset,
        // 0-indexed to start of the probe sequence. Exclusive.
        int endOffset
)
{
    public TargetedRange
    {
        if(!(startOffset >= 0 && startOffset < endOffset))
        {
            throw new IllegalArgumentException("Invalid offsets");
        }
    }

    public static TargetedRange fromStartAndLength(int startOffset, int length)
    {
        return new TargetedRange(startOffset, startOffset + length);
    }

    public static TargetedRange fromRegions(final ChrBaseRegion target, final ChrBaseRegion probeRegion)
    {
        ChrBaseRegion intersection = regionIntersection(target, probeRegion).orElseThrow();
        return fromStartAndLength(intersection.start() - probeRegion.start(), intersection.baseLength());
    }

    public static TargetedRange singleBase(int offset)
    {
        return fromStartAndLength(offset, 1);
    }

    public static TargetedRange wholeRegion(int length)
    {
        return fromStartAndLength(0, length);
    }

    public int baseLength()
    {
        return endOffset - startOffset;
    }
}
