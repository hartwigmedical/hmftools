package com.hartwig.hmftools.svprep;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.reads.PartitionStats;
import com.hartwig.hmftools.svprep.reads.PartitionTask;
import com.hartwig.hmftools.svprep.reads.PartitionThread;
import com.hartwig.hmftools.svprep.reads.ReadFilterType;

public class ChromosomeTask implements AutoCloseable
{
    private final String mChromosome;
    private final SvConfig mConfig;
    private final SpanningReadCache mSpanningReadCache;
    private final ExistingJunctionCache mExistingJunctionCache;
    private final ResultsWriter mWriter;
    private final Queue<PartitionTask> mPartitions;

    private final CombinedStats mCombinedStats;

    public ChromosomeTask(
            final String chromosome, final SvConfig config, final SpanningReadCache spanningReadCache,
            final ExistingJunctionCache existingJunctionCache, final ResultsWriter writer)
    {
        mChromosome = chromosome;
        mConfig = config;
        mSpanningReadCache = spanningReadCache;
        mExistingJunctionCache = existingJunctionCache;
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
    public CombinedStats combinedStats() { return mCombinedStats; }

    public void process()
    {
        int regionCount = mPartitions.size();
        SV_LOGGER.debug("chromosome({}) executing {} regions", mChromosome, regionCount);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(mPartitions.size(), mConfig.Threads); ++i)
        {
            workers.add(new PartitionThread(mChromosome, mConfig, mPartitions, mSpanningReadCache, mExistingJunctionCache, mWriter, mCombinedStats));
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

        SV_LOGGER.info("chromosome({}) {} regions complete, stats: {}",
                mChromosome, regionCount, mCombinedStats.ReadStats.toString());

        if(Arrays.stream(mCombinedStats.ReadStats.ReadFilterCounts).anyMatch(x -> x > 0))
        {
            SV_LOGGER.debug("chromosome({}) filters({})",
                    mChromosome, ReadFilterType.filterCountsToString(mCombinedStats.ReadStats.ReadFilterCounts));
        }

        mSpanningReadCache.logStats();

        if(mConfig.PerfDebug && mCombinedStats.ReadStats.TotalReads > 10000)
        {
            if(SV_LOGGER.isDebugEnabled())
                mCombinedStats.PerfCounters.forEach(x -> x.logIntervalStats(5));
            else
                mCombinedStats.PerfCounters.forEach(x -> x.logStats());
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
