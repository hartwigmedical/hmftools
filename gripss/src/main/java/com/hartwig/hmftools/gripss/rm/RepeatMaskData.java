package com.hartwig.hmftools.gripss.rm;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class RepeatMaskData
{
    public final int Id;
    public final BaseRegion Region;
    public final int SwScore;
    public final char Orientation;
    public final String Repeat;
    public final String ClassType;

    public RepeatMaskData(
            final int id, final BaseRegion region, final int swScore, final char orientation, final String repeat, final String classType)
    {
        Id = id;
        Region = region;
        SwScore = swScore;
        Orientation = orientation;
        Repeat = repeat;
        ClassType = classType;
    }

    public int overlappingBases(final AlignmentData alignment)
    {
        return min(Region.end(), alignment.Region.end()) - max(Region.start(), alignment.Region.start()) + 1;
    }

    public String toString() { return format("%d: region(%s) type(%s: %s)", Id, Region, ClassType, Repeat); }
}
