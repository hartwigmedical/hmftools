package com.hartwig.hmftools.bamtools.fastqbimodalcollapse;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqConfig.registerConfig;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;

import java.util.concurrent.ExecutionException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

// TODO: Multithread
// TODO: Break out into separate classes
// TODO: Reduce reimplementation of common classes.
public class FastqBimodalCollapse
{
    private final FastqBimodalCollapseConfig mConfig;

    public FastqBimodalCollapse(final FastqBimodalCollapseConfig config)
    {
        mConfig = config;
    }

    public void run()
    {
        BT_LOGGER.info("starting FastqBimodalCollapse");

        // TODO:
//        long startTimeMs = System.currentTimeMillis();
//
//        // partition all chromosomes
//        List<ChrBaseRegion> partitions = createPartitions();
//
//        final RemoteReadHandler remoteReadHandler = new RemoteReadHandler(mConfig);
//
//        int numDigits = Integer.toString(mConfig.Threads - 1).length();
//        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + numDigits + "d").build();
//        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
//
//        final ThreadData threadData = new ThreadData(mConfig, remoteReadHandler);
//
//        BT_LOGGER.debug("splitting {} partitions across {} threads", partitions.size(), mConfig.Threads);
//
//        final List<CompletableFuture<Void>> futures = new ArrayList<>();
//
//        addCacheUnmappedReadFutures(futures, remoteReadHandler, executorService);
//
//        // submit each partition to the thread pool
//        for(ChrBaseRegion chrBaseRegion : partitions)
//        {
//            Runnable task = () -> threadData.getPartitionReader().processRegion(chrBaseRegion);
//            futures.add(CompletableFuture.runAsync(task, executorService));
//        }
//
//        // wait for completion
//        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
//
//        BT_LOGGER.info("all partition tasks complete");
//
//        BT_LOGGER.debug("processing cached remote reads");
//        remoteReadHandler.writeRemoteReadsToFastq(executorService, threadData);
//
//        threadData.closePartitionReaders();
//        threadData.closeFastqWriters();
//
//        // might need to combine the fastq files of all the threads
//        mergeThreadFastqFiles(threadData.getAllThreadFastqWriterCaches(), executorService);
//
//        // threading work all done by this point
//        executorService.shutdown();
//
//        // get the total number of reads
//        long totalReads = threadData.getAllThreadFastqWriterCaches().stream().mapToLong(FastqWriterCache::numReadsWritten).sum();
//
//        if(mConfig.PerfDebug)
//        {
//            Iterator<PartitionReader> itr = threadData.getAllThreadPartitionReaders().iterator();
//            PerformanceCounter combinedPerfCounter = itr.next().perfCounter();
//            while(itr.hasNext())
//            {
//                combinedPerfCounter.merge(itr.next().perfCounter());
//            }
//            combinedPerfCounter.logStats();
//        }

        // TODO: Add extra info.
        BT_LOGGER.info("FastqBimodalCollapse complete");
    }



    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        FastqBimodalCollapseConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        FastqBimodalCollapseConfig config = new FastqBimodalCollapseConfig(configBuilder);
        FastqBimodalCollapse fastqBimodalCollapse = new FastqBimodalCollapse(config);
        fastqBimodalCollapse.run();
    }
}
