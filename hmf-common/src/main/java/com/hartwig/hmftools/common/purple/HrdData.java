package com.hartwig.hmftools.common.purple;

public class HrdData
{
    public final int LohSegments;
    public final int SegmentImbalances;
    public final int SegmentBreaks;
    public final HrdStatus Status;

    public static final HrdData INVALID = new HrdData(0, 0, 0, HrdStatus.UNKNOWN);

    public HrdData(final int lohSegments, final int segmentImbalances, final int segmentBreaks, final HrdStatus status)
    {
        LohSegments = lohSegments;
        SegmentImbalances = segmentImbalances;
        SegmentBreaks = segmentBreaks;
        Status = status;
    }
}
