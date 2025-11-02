package com.hartwig.hmftools.redux;

import static java.lang.Math.ceil;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.chromosome.Chromosome.isAltRegionContig;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.formHumanChromosomeRegions;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.perf.TaskQueue;
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

        if(!specificRegions.Regions.isEmpty())
        {
            // only split by thread count if can be done simply
            if(threadCount == specificRegions.Regions.size() && specificRegions.Regions.size() > 1)
            {
                for(int i = 0 ; i < threadCount; ++i)
                {
                    partitionRegions.add(List.of(specificRegions.Regions.get(i)));
                }
            }
            else if(specificRegions.Regions.size() == 1 && threadCount > 1)
            {
                ChrBaseRegion specificRegion = specificRegions.Regions.get(0);
                int intervalLength = (int)ceil(specificRegion.baseLength() / (double)threadCount);
                int regionStart = specificRegion.start();

                for(int i = 0 ; i < threadCount; ++i)
                {
                    int regionEnd = min(regionStart + intervalLength - 1, specificRegion.end());
                    partitionRegions.add(List.of(new ChrBaseRegion(specificRegion.Chromosome, regionStart, regionEnd)));
                    regionStart = regionEnd + 1;
                }
            }
            else
            {
                partitionRegions.addAll(List.of(specificRegions.Regions));
            }

            return partitionRegions;
        }

        List<ChrBaseRegion> genomeRegions = getRefGenomeRegions(specificRegions, refGenomeVersion, refGenome);

        // divvy up standard chromosomes amongst the reads evenly by base length, but for other contigs just add them to the final thread

        // ignore lengths for alt-contigs so they don't impact the partitioning of actual chromosomes and reads
        long totalLength = genomeRegions.stream().mapToLong(x -> isAltRegionContig(x.Chromosome) ? 0 : x.baseLength()).sum();

        long intervalLength = (int)ceil(totalLength / (double)threadCount);

        int chrEndBuffer = (int)round(intervalLength * 0.05);
        long currentLength = 0;
        int nextRegionStart = 1;
        List<ChrBaseRegion> currentRegions = Lists.newArrayList();
        partitionRegions.add(currentRegions);

        for(ChrBaseRegion genomeRegion : genomeRegions)
        {
            if(isAltRegionContig(genomeRegion.Chromosome))
            {
                currentRegions.add(genomeRegion);
                continue;
            }

            nextRegionStart = 1;
            int regionLength = genomeRegion.baseLength();
            int remainingChromosomeLength = regionLength;

            while(currentLength + remainingChromosomeLength >= intervalLength)
            {
                int remainingIntervalLength = (int)(intervalLength - currentLength);
                int regionEnd = nextRegionStart + remainingIntervalLength - 1;

                if(regionLength - regionEnd < chrEndBuffer)
                {
                    regionEnd = regionLength;
                    remainingChromosomeLength = 0;
                }

                currentRegions.add(new ChrBaseRegion(genomeRegion.Chromosome, nextRegionStart, regionEnd));

                currentRegions = Lists.newArrayList();
                partitionRegions.add(currentRegions);
                currentLength = 0;

                if(remainingChromosomeLength == 0)
                    break;

                nextRegionStart = regionEnd + 1;
                remainingChromosomeLength = regionLength - nextRegionStart + 1;
            }

            if(remainingChromosomeLength <= 0)
                continue;

            currentLength += remainingChromosomeLength;

            currentRegions.add(new ChrBaseRegion(genomeRegion.Chromosome, nextRegionStart, regionLength));
        }

        return partitionRegions.stream().filter(x -> !x.isEmpty()).collect(Collectors.toList());
    }

    private static List<ChrBaseRegion> getRefGenomeRegions(
            final SpecificRegions specificRegions, final RefGenomeVersion refGenomeVersion, @Nullable final RefGenomeInterface refGenome)
    {
        if(refGenome == null)
            return formHumanChromosomeRegions(specificRegions, refGenomeVersion);

        CachedRefGenome cachedRefGenome = (CachedRefGenome)refGenome;
        RefGenomeSource refGenomeSource = (RefGenomeSource)cachedRefGenome.refGenome();
        return refGenomeSource.formRefGenomeRegions(specificRegions);
    }
}
