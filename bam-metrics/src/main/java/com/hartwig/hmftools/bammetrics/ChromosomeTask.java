package com.hartwig.hmftools.bammetrics;

import static java.lang.Math.min;

import static com.hartwig.hmftools.bammetrics.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class ChromosomeTask
{
    private final String mChromosome;
    private final BmConfig mConfig;
    private final Queue<PartitionTask> mPartitions;

    private final Metrics mCombinedMetrics;

    public ChromosomeTask(final String chromosome, final BmConfig config)
    {
        mChromosome = chromosome;
        mConfig = config;

        mCombinedMetrics = new Metrics(config.PartitionSize);

        mPartitions = new ConcurrentLinkedQueue<>();
        List<ChrBaseRegion> partitions = partition(chromosome);

        int taskId = 0;
        for(int i = 0; i < partitions.size(); ++i)
        {
            ChrBaseRegion region = partitions.get(i);
            mPartitions.add(new PartitionTask(region, taskId++));
        }
    }

    public Metrics combinedMetrics() { return mCombinedMetrics; }

    public void process()
    {
        int regionCount = mPartitions.size();
        BM_LOGGER.debug("chromosome({}) executing {} regions", mChromosome, regionCount);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(mPartitions.size(), mConfig.Threads); ++i)
        {
            workers.add(new PartitionThread(mChromosome, mConfig, mPartitions, mCombinedMetrics));
        }

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                BM_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        BM_LOGGER.info("chromosome({}) {} regions complete, stats: {}",
                mChromosome, regionCount, mCombinedMetrics.toString());

        if(mConfig.PerfDebug)
        {
            // mCombinedStats.PerfCounters.forEach(x -> x.logStats());
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
}
