package com.hartwig.hmftools.esvee.utils;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class HighDepthRegion extends ChrBaseRegion
{
    public int DepthMin;
    public int DepthMax;
    public int SampleCount;

    public HighDepthRegion(final ChrBaseRegion region)
    {
        super(region.Chromosome, region.start(), region.end());
        DepthMin = 0;
        DepthMax = 0;
        SampleCount = 0;
    }

    public String toString() { return format("region(%s:%d_%d) depth(min=%d max=%d) samples(%d)",
            Chromosome, start(), end(), DepthMin, DepthMax, SampleCount); }
}
