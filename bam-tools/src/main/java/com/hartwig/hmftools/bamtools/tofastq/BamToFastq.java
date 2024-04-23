package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.tofastq.FastqConfig.registerConfig;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

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

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile));

        final SAMFileHeader fileHeader = samReader.getFileHeader();

        for(final SAMSequenceRecord sequenceRecord : fileHeader.getSequenceDictionary().getSequences())
        {
            String chromosome = sequenceRecord.getSequenceName();

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            List<ChrBaseRegion> chrPartitions = Lists.newArrayList();

            if(!mConfig.SpecificChrRegions.Regions.isEmpty())
            {
                chrPartitions.addAll(partitionChromosome(
                        chromosome, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize));
            }
            else
            {
                chrPartitions.addAll(buildPartitions(chromosome, sequenceRecord.getEnd(), mConfig.PartitionSize));
            }

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

        BT_LOGGER.debug("processing partition cache unmatched reads");

        partitionDataStore.processUnmatchedReads(mWriterCache);

        // free up any processing state
        partitionReaders.clear();

        System.gc();

        BT_LOGGER.debug("processing unmapped reads");

        processUnmappedReads();

        BT_LOGGER.debug("writing unpaired reads");

        mWriterCache.writeUnpairedReads();

        mWriterCache.close();

        if(mConfig.PerfDebug)
        {
            double totalPartitionCacheLockTime = partitionDataStore.partitions().stream().mapToDouble(x -> x.totalLockTimeMs()/1000).sum();
            BT_LOGGER.debug("partition store lock time({}}s)", format("%.3f",totalPartitionCacheLockTime));

            PerformanceCounter combinedPerfCounter = partitionReaders.get(0).perfCounter();

            if(partitionReaders.size() > 1)
            {
                for(int i = 1; i < partitionReaders.size(); ++i)
                {
                    combinedPerfCounter.merge(partitionReaders.get(i).perfCounter());
                }
            }

            combinedPerfCounter.logStats();
        }

        BT_LOGGER.info("BamToFastq complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void processUnmappedReads()
    {
        if(mConfig.SpecificChrRegions.hasFilters())
            return;

        UnmappedReads unmappedReads = new UnmappedReads(mConfig, mWriterCache);
        unmappedReads.run();
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
