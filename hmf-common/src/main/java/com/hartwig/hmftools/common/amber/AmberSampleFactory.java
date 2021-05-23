package com.hartwig.hmftools.common.amber;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class AmberSampleFactory
{
    private final int minReadDepth;
    private final NormalHeterozygousFilter heterozygousFilter;

    public AmberSampleFactory(final int minReadDepth, final double minHetAFPercentage, final double maxHetAFPercentage)
    {
        this.minReadDepth = minReadDepth;
        this.heterozygousFilter = new NormalHeterozygousFilter(minHetAFPercentage, maxHetAFPercentage);
    }

    public AmberSample fromBaseDepth(@NotNull final String sample, @NotNull final List<BaseDepth> baseDepths)
    {
        byte[] entries = new byte[baseDepths.size()];
        for(int i = 0; i < baseDepths.size(); i++)
        {
            entries[i] = asByte(baseDepths.get(i));
        }

        return ImmutableAmberSample.builder().sampleId(sample).entries(entries).build();
    }

    public byte asByte(BaseDepth depth)
    {
        if(!depth.isValid() || depth.readDepth() < minReadDepth)
        {
            return AmberSample.DO_NOT_MATCH;
        }

        if(depth.refSupport() == depth.readDepth())
        {
            return (byte) 1;
        }

        if(heterozygousFilter.test(depth))
        {
            return (byte) 2;
        }

        if(depth.altSupport() == depth.readDepth())
        {
            return (byte) 3;
        }

        return AmberSample.DO_NOT_MATCH;
    }

}
