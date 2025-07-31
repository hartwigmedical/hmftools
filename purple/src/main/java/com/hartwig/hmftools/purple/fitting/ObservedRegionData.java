package com.hartwig.hmftools.purple.fitting;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.region.FittingRegion;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public class ObservedRegionData
{
    public final FittingRegion Region;
    public final List<SomaticVariant> Variants;

    public ObservedRegionData(final FittingRegion region)
    {
        Region = region;
        Variants = Lists.newArrayListWithExpectedSize(2);
    }

    public void addVariant(final SomaticVariant variant) { Variants.add(variant); }
}
