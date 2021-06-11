package com.hartwig.hmftools.sage.pipeline;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
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
    public List<GenomeRegion> partition(String contig)
    {
        return partition(contig, 1, mRefGenome.getSequence(contig).length());
    }

    @NotNull
    public List<GenomeRegion> partition(String contig, int minPosition, int maxPosition)
    {
        final List<GenomeRegion> results = Lists.newArrayList();

        int dynamicSliceSize = maxPosition / Math.min(mConfig.Threads, 4) + 1;

        final int regionSliceSize = Math.min(dynamicSliceSize, mConfig.RegionSliceSize);
        for(int i = 0; ; i++)
        {
            int start = minPosition + i * regionSliceSize;
            int end = Math.min(start + regionSliceSize - 1, maxPosition);
            results.add(GenomeRegions.create(contig, start, end));
            if(end >= maxPosition)
            {
                break;
            }
        }
        return results;
    }

}
