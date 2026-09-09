package com.hartwig.hmftools.common.region;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class HighDepthRegion extends ChrBaseRegion
{
    public int DepthMin;
    public int DepthMax;
    public int DepthAvg;
    public int SampleCount;

    public HighDepthRegion(final ChrBaseRegion region)
    {
        super(region.Chromosome, region.start(), region.end());
        DepthMin = 0;
        DepthMax = 0;
        DepthAvg = 0;
        SampleCount = 0;
    }

    public HighDepthRegion(
            final String chromosome, final int posStart, final int posEnd,
            final int depthMin, final int depthMax, final int depthAvg, final int sampleCount)
    {
        super(chromosome, posStart, posEnd);
        DepthMin = depthMin;
        DepthMax = depthMax;
        DepthAvg = depthAvg;
        SampleCount = sampleCount;
    }

    // file loading
    public static final String FLD_SAMPLE_COUNT = "SampleCount";
    public static final String FLD_DEPTH_MIN = "DepthMin";
    public static final String FLD_DEPTH_MAX = "DepthMax";
    public static final String FLD_DEPTH_AVG = "DepthAvg";

    public String toString() { return format("region(%s:%d_%d) depth(min=%d max=%d avg=%d) samples(%d)",
            Chromosome, start(), end(), DepthMin, DepthMax, DepthAvg, SampleCount); }
}
