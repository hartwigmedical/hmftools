package com.hartwig.hmftools.common.mappability;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.BaseRegion;

public class RegionQuality extends BaseRegion
{
    public final double Quality;

    public RegionQuality(final BaseRegion region, double quality)
    {
        super(region.start(), region.end());
        Quality = quality;
    }

    public String toString() { return format("%s: %.3f", super.toString(), Quality); }
}
