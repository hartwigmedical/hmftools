package com.hartwig.hmftools.common.purple;

public class HrdData
{
    public final int LohSegments;
    public final int SegmentImbalances;
    public final int SegmentBreaks;
    public final HrdStatus Status;

    public HrdData(final int lohSegments, final int segmentImbalances, final int segmentBreaks, final HrdStatus status)
    {
        LohSegments = lohSegments;
        SegmentImbalances = segmentImbalances;
        SegmentBreaks = segmentBreaks;
        Status = status;
    }
}
