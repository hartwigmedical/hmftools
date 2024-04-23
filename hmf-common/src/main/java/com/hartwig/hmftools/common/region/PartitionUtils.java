package com.hartwig.hmftools.common.region;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;

import static org.apache.commons.lang3.ObjectUtils.max;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import htsjdk.samtools.SAMSequenceRecord;

public class PartitionUtils
{
    public static List<ChrBaseRegion> partitionChromosome(
            final String chromosome, final RefGenomeVersion refGenomeVersion, final List<ChrBaseRegion> specificRegions, int partitionSize)
    {
        if(!specificRegions.isEmpty())
        {
            final List<ChrBaseRegion> chromosomeRegions = specificRegions.stream()
                    .filter(x -> x.Chromosome.equals(chromosome)).collect(Collectors.toList());

            return partitionSpecificRegions(chromosomeRegions, partitionSize);
        }
        else
        {
            return partitionChromosome(chromosome, refGenomeVersion, partitionSize);
        }
    }

    public static List<ChrBaseRegion> partitionChromosome(final String chromosome, final RefGenomeVersion refGenomeVersion, int partitionSize)
    {
        RefGenomeCoordinates refGenomeCoords = refGenomeVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        int chromosomeLength = refGenomeCoords.length(stripChrPrefix(chromosome));
        return buildPartitions(refGenomeVersion.versionedChromosome(chromosome), chromosomeLength, partitionSize);
    }

    public static List<ChrBaseRegion> partitionChromosome(final SAMSequenceRecord refSequence, int partitionSize)
    {
        int chromosomeLength = refSequence.getSequenceLength();
        return buildPartitions(refSequence.getContig(), chromosomeLength, partitionSize);
    }

    public static List<ChrBaseRegion> partitionSpecificRegions(final List<ChrBaseRegion> specificRegions, int partitionSize)
    {
        List<ChrBaseRegion> partitions = Lists.newArrayList();

        for(ChrBaseRegion region : specificRegions)
        {
            // break regions across partition lines if required
            int partitionIndexStart = region.start() / partitionSize;
            int partitionIndexEnd = region.end() / partitionSize;

            if(partitionIndexStart == partitionIndexEnd)
            {
                partitions.add(region);
            }
            else
            {
                for(int i = partitionIndexStart; i <= partitionIndexEnd; ++i)
                {
                    int partitionStart = i * partitionSize;
                    int partitionEnd = partitionStart + partitionSize - 1;

                    partitions.add(new ChrBaseRegion(
                            region.Chromosome, max(max(partitionStart, region.start()), 1), min(partitionEnd, region.end())));
                }
            }
        }

        return partitions;
    }

    public static List<ChrBaseRegion> buildPartitions(final String chromosome, int maxPosition, int partitionSize)
    {
        final List<ChrBaseRegion> partitions = Lists.newArrayList();

        int minPosition = 0; // actually starts at 1 for first partition

        for(int i = 0; ; i++)
        {
            int start = minPosition + i * partitionSize;
            int end = min(start + partitionSize - 1, maxPosition);
            partitions.add(new ChrBaseRegion(chromosome, max(start, 1), end));

            if(end >= maxPosition)
                break;
        }

        return partitions;
    }
}
