package com.hartwig.hmftools.purple.hrd;

public class HrdData
{
    public final int LohSegments;
    public final int SegmentImbalances;
    public final int SegmentBreaks;

    public HrdData(final int lohSegments, final int segmentImbalances, final int segmentBreaks)
    {
        LohSegments = lohSegments;
        SegmentImbalances = segmentImbalances;
        SegmentBreaks = segmentBreaks;
    }
}
