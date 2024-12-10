package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import com.hartwig.hmftools.common.genome.region.Orientation;

public class SupplementaryReadInfo
{
    public final int UnclippedPosition;
    public final Orientation Orient;

    public SupplementaryReadInfo(final int unclippedPosition, final Orientation orient)
    {
        UnclippedPosition = unclippedPosition;
        Orient = orient;
    }

    public String toString()
    {
        return format("%d:%d", UnclippedPosition, Orient.asByte());
    }
}
