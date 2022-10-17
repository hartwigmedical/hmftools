package com.hartwig.hmftools.svprep.tools;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.compress.utils.Lists;

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
