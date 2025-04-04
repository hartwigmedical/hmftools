package com.hartwig.hmftools.bamtools.common;

import static java.lang.Math.ceil;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.checker.CheckConfig;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionTask
{
    public final ChrBaseRegion Region;
    public final int TaskId;

    public PartitionTask(final ChrBaseRegion region, final int taskId)
    {
        Region = region;
        TaskId = taskId;
    }

    public static List<ChrBaseRegion> splitRegionsIntoPartitions(
            final String bamFile, final String refGenomeFile, final int threads, final SpecificRegions specificRegions, final int partitionSize)
    {
        List<ChrBaseRegion> partitionRegions = Lists.newArrayList();

        if(!specificRegions.Regions.isEmpty())
        {
            // only split by thread count if can be done simply
            if(specificRegions.Regions.size() == 1)
            {
                if(threads > 1)
                {
                    ChrBaseRegion specificRegion = specificRegions.Regions.get(0);

                    int regionCount = (int)ceil(specificRegion.baseLength() / (double)partitionSize);

                    // int intervalLength = (int)ceil(specificRegion.baseLength() / (double)threads);
                    int regionStart = specificRegion.start();

                    for(int i = 0 ; i < regionCount; ++i)
                    {
                        int regionEnd = min(regionStart + partitionSize - 1, specificRegion.end());
                        partitionRegions.add(new ChrBaseRegion(specificRegion.Chromosome, regionStart, regionEnd));
                        regionStart = regionEnd + 1;
                    }
                }
                else
                {
                    partitionRegions.add(specificRegions.Regions.get(0));
                }
            }
            else
            {
                for(ChrBaseRegion region : specificRegions.Regions)
                {
                    int regionCount = (int)ceil(region.baseLength() / (double)partitionSize);

                    // int intervalLength = (int)ceil(specificRegion.baseLength() / (double)threads);
                    int regionStart = region.start();

                    for(int i = 0 ; i < regionCount; ++i)
                    {
                        int regionEnd = min(regionStart + partitionSize - 1, region.end());
                        partitionRegions.add(new ChrBaseRegion(region.Chromosome, regionStart, regionEnd));
                        regionStart = regionEnd + 1;
                    }
                }
            }

            return partitionRegions;
        }

        SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(refGenomeFile))
                .open(new File(bamFile));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        for(final SAMSequenceRecord sequenceRecord : fileHeader.getSequenceDictionary().getSequences())
        {
            String chromosome = sequenceRecord.getSequenceName();

            if(specificRegions != null && specificRegions.excludeChromosome(chromosome))
                continue;

            partitionRegions.addAll(buildPartitions(chromosome, sequenceRecord.getEnd(), partitionSize));
        }

        return partitionRegions;
    }

}
