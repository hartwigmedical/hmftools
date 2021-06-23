package com.hartwig.hmftools.sage.pipeline;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class ChromosomePartition
{
    private final SageConfig mConfig;
    private final ReferenceSequenceFile mRefGenome;

    public ChromosomePartition(final SageConfig config, final ReferenceSequenceFile refGenome)
    {
        mConfig = config;
        mRefGenome = refGenome;
    }

    @NotNull
    public List<BaseRegion> partition(String contig)
    {
        return partition(contig, 1, mRefGenome.getSequence(contig).length());
    }

    @NotNull
    public List<BaseRegion> partition(String contig, int minPosition, int maxPosition)
    {
        final List<BaseRegion> results = Lists.newArrayList();

        int dynamicSliceSize = maxPosition / Math.min(mConfig.Threads, 4) + 1;
        int regionSliceSize = Math.min(dynamicSliceSize, mConfig.RegionSliceSize);

        for(int i = 0; ; i++)
        {
            int start = minPosition + i * regionSliceSize;
            int end = Math.min(start + regionSliceSize - 1, maxPosition);
            results.add(new BaseRegion(contig, start, end));

            if(end >= maxPosition)
                break;
        }
        return results;
    }

}
