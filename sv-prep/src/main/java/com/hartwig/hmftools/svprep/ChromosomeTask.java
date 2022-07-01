package com.hartwig.hmftools.svprep;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class ChromosomeTask implements AutoCloseable
{
    private final String mChromosome;
    private final SvConfig mConfig;
    private final ResultsWriter mWriter;
    private final Queue<PartitionTask> mPartitions;

    private final CombinedStats mCombinedStats;

    public ChromosomeTask(final String chromosome, final SvConfig config, final ResultsWriter writer)
    {
        mChromosome = chromosome;
        mConfig = config;
        mWriter = writer;

        mCombinedStats = new CombinedStats();

        mPartitions = new ConcurrentLinkedQueue<>();
        List<ChrBaseRegion> partitions = partition(chromosome);

        int taskId = 0;
        for(int i = 0; i < partitions.size(); ++i)
        {
            ChrBaseRegion region = partitions.get(i);
            mPartitions.add(new PartitionTask(region, taskId++));
        }
    }

    public String chromosome()
    {
        return mChromosome;
    }

    public void process()
    {
        int regionCount = mPartitions.size();
        SV_LOGGER.info("chromosome({}) executing {} regions", mChromosome, regionCount);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(mPartitions.size(), mConfig.Threads); ++i)
        {
            workers.add(new PartitionThread(mChromosome, mConfig, mPartitions, mWriter, mCombinedStats));
        }

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                SV_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        SV_LOGGER.debug("chromosome({}) {} regions complete, stats: {}",
                mChromosome, regionCount, mCombinedStats.ReadStats.toString());

        /*
        if(mConfig.logPerfStats())
        {
            mRegionResults.logPerfCounters();
            SV_LOGGER.debug("chromosome({}) max memory({})", mChromosome, mRegionResults.maxMemoryUsage());
        }
         */

        SV_LOGGER.info("chromosome({}) analysis complete", mChromosome);

        if(SV_LOGGER.isDebugEnabled())
        {
            mCombinedStats.PerfCounter.logStats();
        }
    }

    private List<ChrBaseRegion> partition(final String chromosome)
    {
        if(!mConfig.SpecificRegions.isEmpty())
        {
            List<ChrBaseRegion> partitions = Lists.newArrayList();

            for(ChrBaseRegion region : mConfig.SpecificRegions)
            {
                if(region.Chromosome.equals(chromosome))
                {
                    partitions.addAll(partition(chromosome, region.start() ,region.end()));
                }
            }

            return partitions;
        }

        RefGenomeCoordinates refGenomeCoords = mConfig.RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        int chromosomeLength = refGenomeCoords.length(stripChrPrefix(chromosome));
        return partition(chromosome, 1, chromosomeLength);
    }

    private List<ChrBaseRegion> partition(final String chromosome, int minPosition, int maxPosition)
    {
        final List<ChrBaseRegion> partitions = Lists.newArrayList();

        for(int i = 0; ; i++)
        {
            int start = minPosition + i * mConfig.PartitionSize;
            int end = min(start + mConfig.PartitionSize - 1, maxPosition);
            partitions.add(new ChrBaseRegion(chromosome, start, end));

            if(end >= maxPosition)
                break;
        }

        return partitions;
    }

    /*
    public int maxMemoryUsage() { return mRegionResults.maxMemoryUsage(); }
    */

    @Override
    public void close()
    {
        // mRefGenome.close();
    }
}
