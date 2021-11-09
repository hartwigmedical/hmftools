package com.hartwig.hmftools.purple.drivers;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class DeletionRegionFrequency
{
    public final BaseRegion Region;
    public int Frequency;

    public DeletionRegionFrequency(final BaseRegion region, final int frequency)
    {
        Region = region;
        Frequency = frequency;
    }
}
