package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class Partition extends ChrBaseRegion
{
    public final List<TargetRegionData> TargetRegions;

    public Partition(final ChrBaseRegion region)
    {
        super(region.Chromosome, region.start(), region.end());
        TargetRegions = Lists.newArrayList();
    }

    public String toString() { return format("region(%s) targeted(%d)", super.toString(), TargetRegions.size()); }
}
