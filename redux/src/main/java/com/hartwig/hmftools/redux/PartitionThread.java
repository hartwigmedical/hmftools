package com.hartwig.hmftools.redux;

import static java.lang.Math.ceil;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.refGenomeCoordinates;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.redux.unmap.TaskQueue;
import com.hartwig.hmftools.redux.write.FileWriterCache;
import com.hartwig.hmftools.redux.write.PartitionInfo;

import org.jetbrains.annotations.Nullable;

public class PartitionThread extends Thread
{
    private final FileWriterCache mFileWriterCache;
    private final BamReader mBamReader;
    private final TaskQueue mPartitions;

    private final PartitionReader mPartitionReader;

    public PartitionThread(
            final ReduxConfig config, final List<String> inputBams, final TaskQueue partitions, final FileWriterCache fileWriterCache)
    {
        mFileWriterCache = fileWriterCache;
        mPartitions = partitions;

        mBamReader = new BamReader(inputBams, config.RefGenomeFile);

        mPartitionReader = new PartitionReader(config, mBamReader);
    }

    public PartitionReader partitionReader() { return mPartitionReader; }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionInfo partition = (PartitionInfo)mPartitions.removeItem();

                mPartitionReader.processPartition(partition);

                mFileWriterCache.addCompletedPartition(partition);
            }
            catch(NoSuchElementException e)
            {
                RD_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public static List<List<ChrBaseRegion>> splitRegionsIntoPartitions(
            final SpecificRegions specificRegions, final int threadCount, final RefGenomeVersion refGenomeVersion,
            @Nullable final RefGenomeInterface refGenome)
    {
        List<List<ChrBaseRegion>> partitionRegions = Lists.newArrayList();

        List<ChrBaseRegion> inputRegions = Lists.newArrayList();

        if(!specificRegions.Regions.isEmpty())
        {
            inputRegions.addAll(specificRegions.Regions);
        }
        else
        {
            inputRegions.addAll(getRefGenomeRegions(specificRegions, refGenomeVersion, refGenome));
        }

        long totalLength = inputRegions.stream().mapToLong(x -> x.baseLength()).sum();
        long intervalLength = (int)ceil(totalLength / (double)threadCount);

        int chrEndBuffer = (int)round(intervalLength * 0.05);
        long currentLength = 0;
        int nextRegionStart = 1;
        List<ChrBaseRegion> currentRegions = Lists.newArrayList();
        partitionRegions.add(currentRegions);

        for(ChrBaseRegion inputRegion : inputRegions)
        {
            nextRegionStart = 1;
            int chromosomeLength = inputRegion.baseLength();
            int remainingChromosomeLength = chromosomeLength;

            while(currentLength + remainingChromosomeLength >= intervalLength)
            {
                int remainingIntervalLength = (int)(intervalLength - currentLength);
                int regionEnd = nextRegionStart + remainingIntervalLength - 1;

                if(chromosomeLength - regionEnd < chrEndBuffer)
                {
                    regionEnd = chromosomeLength;
                    remainingChromosomeLength = 0;
                }

                currentRegions.add(new ChrBaseRegion(inputRegion.Chromosome, nextRegionStart, regionEnd));

                currentRegions = Lists.newArrayList();
                partitionRegions.add(currentRegions);
                currentLength = 0;

                if(remainingChromosomeLength == 0)
                    break;

                nextRegionStart = regionEnd + 1;
                remainingChromosomeLength = chromosomeLength - nextRegionStart + 1;
            }

            if(remainingChromosomeLength <= 0)
                continue;

            currentLength += remainingChromosomeLength;

            currentRegions.add(new ChrBaseRegion(inputRegion.Chromosome, nextRegionStart, chromosomeLength));
        }

        return partitionRegions.stream().filter(x -> !x.isEmpty()).collect(Collectors.toList());
    }

    private static List<ChrBaseRegion> getRefGenomeRegions(
            final SpecificRegions specificRegions, final RefGenomeVersion refGenomeVersion, @Nullable final RefGenomeInterface refGenome)
    {
        if(refGenome == null)
            return humanChromosomeRegions(specificRegions, refGenomeVersion);

        List<ChrBaseRegion> inputRegions = Lists.newArrayList();

        for(Map.Entry<String,Integer> contigEntry : refGenome.chromosomeLengths().entrySet())
        {
            String contig = contigEntry.getKey();

            if(specificRegions.excludeChromosome(contig))
                continue;

            int contigLength = contigEntry.getValue();

            inputRegions.add(new ChrBaseRegion(contig, 1, contigLength));
        }

        Collections.sort(inputRegions);

        return inputRegions;
    }

    private static List<ChrBaseRegion> humanChromosomeRegions(final SpecificRegions specificRegions, final RefGenomeVersion refGenomeVersion)
    {
        List<ChrBaseRegion> inputRegions = Lists.newArrayList();

        RefGenomeCoordinates refGenomeCoordinates = refGenomeCoordinates(refGenomeVersion);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = refGenomeVersion.versionedChromosome(chromosome.toString());

            if(specificRegions.excludeChromosome(chromosomeStr))
                continue;

            inputRegions.add(new ChrBaseRegion(chromosomeStr, 1, refGenomeCoordinates.Lengths.get(chromosome)));
        }

        return inputRegions;
    }

}
