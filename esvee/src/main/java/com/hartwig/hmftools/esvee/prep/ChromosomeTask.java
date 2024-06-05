package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.CombinedStats;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterType;

public class ChromosomeTask
{
    private final String mChromosome;
    private final PrepConfig mConfig;
    private final SpanningReadCache mSpanningReadCache;
    private final ResultsWriter mWriter;
    private final Queue<ChrBaseRegion> mPartitions;

    private final CombinedStats mCombinedStats;

    public ChromosomeTask(
            final String chromosome, final PrepConfig config, final SpanningReadCache spanningReadCache, final ResultsWriter writer)
    {
        mChromosome = chromosome;
        mConfig = config;
        mSpanningReadCache = spanningReadCache;
        mWriter = writer;

        mCombinedStats = new CombinedStats();

        mPartitions = new ConcurrentLinkedQueue<>();
        List<ChrBaseRegion> partitions = partition(chromosome);

        int taskId = 0;
        for(int i = 0; i < partitions.size(); ++i)
        {
            ChrBaseRegion region = partitions.get(i);
            mPartitions.add(region);
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
            workers.add(new PartitionThread(mChromosome, mConfig, mPartitions, mSpanningReadCache, mWriter, mCombinedStats));
        }

        if(!runThreadTasks(workers))
            System.exit(1);

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
        if(!mConfig.SpecificChrRegions.Regions.isEmpty())
        {
            List<ChrBaseRegion> partitions = Lists.newArrayList();

            for(ChrBaseRegion region : mConfig.SpecificChrRegions.Regions)
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
}
