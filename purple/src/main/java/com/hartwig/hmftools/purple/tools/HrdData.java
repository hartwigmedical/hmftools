package com.hartwig.hmftools.purple.tools;

public class HrdData
{
    public final int LohSegments;
    public final int SegmentImbalances;
    public final double SegmentBreaks;

    public HrdData(final int lohSegments, final int segmentImbalances, final double segmentBreaks)
    {
        LohSegments = lohSegments;
        SegmentImbalances = segmentImbalances;
        SegmentBreaks = segmentBreaks;
    }

    public double score() { return LohSegments + SegmentImbalances + SegmentBreaks; }
}
