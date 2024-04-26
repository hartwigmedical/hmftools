package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.tofastq.FastqConfig.CHR_UNMAPPED;
import static com.hartwig.hmftools.bamtools.tofastq.FastqConfig.registerConfig;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
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

    public void run() throws ExecutionException, InterruptedException
    {
        BT_LOGGER.info("starting BamToFastq");

        long startTimeMs = System.currentTimeMillis();

        // partition all chromosomes
        List<ChrBaseRegion> partitions = createPartitions();

        final RemoteReadHandler remoteReadHandler = new RemoteReadHandler(mConfig, mWriterCache);

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%02d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        final ThreadData threadData = new ThreadData(mConfig, mWriterCache, remoteReadHandler);
        int partitionCount = partitions.size();

        BT_LOGGER.debug("splitting {} partitions across {} threads", partitionCount, mConfig.Threads);

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        if(!mConfig.SpecificChrRegions.hasFilters() || mConfig.SpecificChrRegions.Chromosomes.contains(CHR_UNMAPPED))
        {
            // write all unmapped reads to hash bams
            futures.add(CompletableFuture.runAsync(remoteReadHandler::cacheAllUnmappedReads, executorService));
        }

        // submit each partition to the thread pool
        for(ChrBaseRegion chrBaseRegion : partitions)
        {
            Runnable task = () -> threadData.getPartitionReader().processRegion(chrBaseRegion);
            futures.add(CompletableFuture.runAsync(task, executorService));
        }

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();

        BT_LOGGER.info("all partition tasks complete");

        System.gc();

        BT_LOGGER.debug("processing cached unmatched reads");
        remoteReadHandler.writeRemoteReadsToFastq(executorService, threadData);

        // threading work all done by this point
        executorService.shutdown();
        threadData.close();

        BT_LOGGER.debug("writing unpaired reads");
        mWriterCache.writeUnpairedReads();
        mWriterCache.close();

        if(mConfig.PerfDebug)
        {
            PerformanceCounter combinedPerfCounter = threadData.getPartitionReaders().get(0).perfCounter();

            if(threadData.getPartitionReaders().size() > 1)
            {
                for(int i = 1; i < threadData.getPartitionReaders().size(); ++i)
                {
                    combinedPerfCounter.merge(threadData.getPartitionReaders().get(i).perfCounter());
                }
            }

            combinedPerfCounter.logStats();
        }

        BT_LOGGER.info("BamToFastq complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<ChrBaseRegion> createPartitions()
    {
        final SAMFileHeader fileHeader;
        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile)))
        {
            fileHeader = samReader.getFileHeader();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        List<ChrBaseRegion> partitions = new ArrayList<>();

        for(final SAMSequenceRecord sequenceRecord : fileHeader.getSequenceDictionary().getSequences())
        {
            String chromosome = sequenceRecord.getSequenceName();

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            if(!mConfig.SpecificChrRegions.Regions.isEmpty())
            {
                partitions.addAll(partitionChromosome(
                        chromosome, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize));
            }
            else
            {
                partitions.addAll(buildPartitions(chromosome, sequenceRecord.getEnd(), mConfig.PartitionSize));
            }
        }

        return partitions;
    }

    public static void main(final String[] args) throws ExecutionException, InterruptedException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        // set all thread exception handler
        // we must do this otherwise unhandled exception in other threads will cause the program
        // to hang and never return. This will cause the VM to hang around forever.
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            BT_LOGGER.error("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        BamToFastq bamToFastq = new BamToFastq(configBuilder);
        bamToFastq.run();
    }
}
