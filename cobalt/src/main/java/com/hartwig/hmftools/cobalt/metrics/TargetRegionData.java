package com.hartwig.hmftools.cobalt.metrics;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class TargetRegionData extends ChrBaseRegion
{
    public final FragmentGcMap FragmentGcCounts;

    public TargetRegionData(final String chromosome, final BaseRegion region)
    {
        super(chromosome, region.start(), region.end());

        FragmentGcCounts = new FragmentGcMap();
    }

    public String toString() { return super.toString(); }
}
