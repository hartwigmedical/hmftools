package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.config.SageConfig;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public class ChromosomePartition
{
    private final SageConfig mConfig;
    private final ReferenceSequenceFile mRefGenome;
    private final RefGenomeCoordinates mRefCoords;

    public ChromosomePartition(final SageConfig config, final ReferenceSequenceFile refGenome, final RefGenomeCoordinates refCoords)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mRefCoords = refCoords;
    }

    public List<ChrBaseRegion> partition(final String chromosome)
    {
        if(!mConfig.SpecificRegions.isEmpty())
        {
            List<ChrBaseRegion> chrRegions = mConfig.SpecificRegions.stream().filter(x -> x.Chromosome.equals(chromosome)).collect(Collectors.toList());
            return partitionRegions(chrRegions);
        }

        int chromosomeLength = mRefCoords.lengths().get(chromosome).intValue();
        return partition(chromosome, 1, chromosomeLength);
    }

    public List<ChrBaseRegion> partition(final String chromosome, int minPosition, int maxPosition)
    {
        final List<ChrBaseRegion> partitions = Lists.newArrayList();

        int dynamicSliceSize = maxPosition / Math.min(mConfig.Threads, 4) + 1;
        int regionSliceSize = Math.min(dynamicSliceSize, mConfig.RegionSliceSize);

        for(int i = 0; ; i++)
        {
            int start = minPosition + i * regionSliceSize;
            int end = Math.min(start + regionSliceSize - 1, maxPosition);
            partitions.add(new ChrBaseRegion(chromosome, start, end));

            if(end >= maxPosition)
                break;
        }
        return partitions;
    }

    public List<ChrBaseRegion> partitionRegions(final List<ChrBaseRegion> regions)
    {
        List<ChrBaseRegion> partitions = Lists.newArrayList();

        int totalRegionLength = regions.stream().mapToInt(x -> x.baseLength()).sum();
        int dynamicSliceSize = totalRegionLength / Math.min(mConfig.Threads, 4) + 1;

        for(ChrBaseRegion region : regions)
        {
            if(region.baseLength() < mConfig.RegionSliceSize)
            {
                partitions.add(region);
            }
            else
            {
                int minPosition = region.start();
                int maxPosition = region.end();
                int regionSliceSize = Math.min(dynamicSliceSize, mConfig.RegionSliceSize);

                for(int i = 0; ; i++)
                {
                    int start = minPosition + i * regionSliceSize;

                    int end = start + regionSliceSize - 1;

                    if(end > maxPosition)
                        end = maxPosition;
                    else if(maxPosition - end < 0.5 * regionSliceSize) // extend the last segment if within 50% of the end
                        end = maxPosition;

                    partitions.add(new ChrBaseRegion(region.Chromosome, start, end));

                    if(end >= maxPosition)
                        break;
                }
            }
        }

        return partitions;
    }

}
