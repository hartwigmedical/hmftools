package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.tofastq.FastqConfig.registerConfig;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMRecord;

public class BamToFastq
{
    private final FastqConfig mConfig;
    private final FastqWriterCache mWriterCache;

    public BamToFastq(final ConfigBuilder configBuilder)
    {
        mConfig = new FastqConfig(configBuilder);
        mWriterCache = new FastqWriterCache(mConfig);
    }

    public void run()
    {
        BT_LOGGER.info("starting BamToFastq");

        long startTimeMs = System.currentTimeMillis();

        PartitionDataStore partitionDataStore = new PartitionDataStore(mConfig);

        // partition all chromosomes
        Queue<ChrBaseRegion> partitions = new ConcurrentLinkedQueue<>();

        // TO-DO: loop through BAM header instead of human chromosomes

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            List<ChrBaseRegion> chrPartitions = partitionChromosome(
                    chromosomeStr, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize);

            partitions.addAll(chrPartitions);
        }

        List<PartitionReader> partitionReaders = Lists.newArrayList();
        List<Thread> workers = new ArrayList<>();

        int partitionCount = partitions.size();

        for(int i = 0; i < min(partitionCount, mConfig.Threads); ++i)
        {
            PartitionReader partitionThread = new PartitionReader(mConfig, partitions, mWriterCache, partitionDataStore);
            partitionReaders.add(partitionThread);
            workers.add(partitionThread);
        }

        BT_LOGGER.debug("splitting {} partitions across {} threads", partitionCount, partitionReaders.size());

        if(!runThreadTasks(workers))
            System.exit(1);

        BT_LOGGER.info("all partition tasks complete");

        List<SAMRecord> unmatchedReads = Lists.newArrayList();

        for(PartitionData partitionData : partitionDataStore.partitions())
        {
            unmatchedReads.addAll(partitionData.unmatchedReads());
        }

        if(!unmatchedReads.isEmpty())
        {
            // TO-DO: what to do with these? suggest a bug somewhere

            BT_LOGGER.info("unmatched read count {}", unmatchedReads.size());

            for(int i = 0; i < unmatchedReads.size(); ++i)
            {
                SAMRecord read = unmatchedReads.get(i);

                BT_LOGGER.warn("unmatched read: {}", readToString(read));
            }
        }

        PerformanceCounter combinedPerfCounter = partitionReaders.get(0).perfCounter();

        if(partitionReaders.size() > 1)
        {
            for(int i = 1; i < partitionReaders.size(); ++i)
            {
                combinedPerfCounter.merge(partitionReaders.get(i).perfCounter());
            }
        }

        combinedPerfCounter.logStats();

        // free up any processing state
        partitionReaders.clear();

        System.gc();

        processUnmappedReads();

        mWriterCache.close();

        BT_LOGGER.info("BamToFastq complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void processUnmappedReads()
    {
        if(mConfig.SpecificChrRegions.hasFilters())
            return;

        UnmappedReads unmappedReads = new UnmappedReads(mConfig, mWriterCache);
        unmappedReads.run();

        BT_LOGGER.info("processed {} unmapped reads", unmappedReads.unmappedReadCount());
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BamToFastq bamToFastq = new BamToFastq(configBuilder);
        bamToFastq.run();
    }
}
