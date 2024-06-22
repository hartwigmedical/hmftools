package com.hartwig.hmftools.redux.merge.highdepth;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class HighDepthRegion extends ChrBaseRegion
{
    public int DepthMin;
    public int DepthMax;
    public long BaseVolume;

    public HighDepthRegion(final ChrBaseRegion region)
    {
        super(region.Chromosome, region.start(), region.end());
        DepthMin = 0;
        DepthMax = 0;
        BaseVolume = 0;
    }

    public void merge(final HighDepthRegion other)
    {
        DepthMin = Math.min(DepthMin, other.DepthMin);
        DepthMax = Math.max(DepthMax, other.DepthMax);
        BaseVolume += other.BaseVolume;
        setStart(Math.min(start(), other.start()));
        setEnd(Math.max(end(), other.end()));
    }

    public float depthAvg()
    {
        return 1.0f * BaseVolume / baseLength();
    }

    public String toString()
    {
        return format("region(%s:%d_%d) depth(min=%d max=%d avg=%.2f)",
                Chromosome, start(), end(), DepthMin, DepthMax, depthAvg());
    }
}
