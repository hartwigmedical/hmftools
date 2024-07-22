package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqConfig.CHR_UNMAPPED;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqConfig.registerConfig;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.R1;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.R2;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.UNPAIRED;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.formFilename;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.mergeFastqs;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.Level;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;

public class BamToFastq
{
    private final ToFastqConfig mConfig;

    public BamToFastq(final ConfigBuilder configBuilder)
    {
        mConfig = new ToFastqConfig(configBuilder);
    }

    public void run() throws ExecutionException, InterruptedException
    {
        BT_LOGGER.info("starting BamToFastq");

        long startTimeMs = System.currentTimeMillis();

        // partition all chromosomes
        List<ChrBaseRegion> partitions = createPartitions();

        final RemoteReadHandler remoteReadHandler = new RemoteReadHandler(mConfig);

        int numDigits = Integer.toString(mConfig.Threads - 1).length();
        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + numDigits + "d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        final ThreadData threadData = new ThreadData(mConfig, remoteReadHandler);

        BT_LOGGER.debug("splitting {} partitions across {} threads", partitions.size(), mConfig.Threads);

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        addCacheUnmappedReadFutures(futures, remoteReadHandler, executorService);

        // submit each partition to the thread pool
        for(ChrBaseRegion chrBaseRegion : partitions)
        {
            Runnable task = () -> threadData.getPartitionReader().processRegion(chrBaseRegion);
            futures.add(CompletableFuture.runAsync(task, executorService));
        }

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();

        BT_LOGGER.info("all partition tasks complete");

        BT_LOGGER.debug("processing cached remote reads");
        remoteReadHandler.writeRemoteReadsToFastq(executorService, threadData);

        threadData.closePartitionReaders();
        threadData.closeFastqWriters();

        // might need to combine the fastq files of all the threads
        mergeThreadFastqFiles(threadData.getAllThreadFastqWriterCaches(), executorService);

        // threading work all done by this point
        executorService.shutdown();

        // get the total number of reads
        long totalReads = threadData.getAllThreadFastqWriterCaches().stream().mapToLong(FastqWriterCache::numReadsWritten).sum();

        if(mConfig.PerfDebug)
        {
            Iterator<PartitionReader> itr = threadData.getAllThreadPartitionReaders().iterator();
            PerformanceCounter combinedPerfCounter = itr.next().perfCounter();
            while(itr.hasNext())
            {
                combinedPerfCounter.merge(itr.next().perfCounter());
            }
            combinedPerfCounter.logStats();
        }

        BT_LOGGER.printf(Level.INFO, "BamToFastq complete, total reads(%,d), time taken(%s mins)",
                totalReads, runTimeMinsStr(startTimeMs));
    }

    private List<ChrBaseRegion> createPartitions()
    {
        final SAMFileHeader fileHeader;
        try(SamReader samReader = ToFastqUtils.openSamReader(mConfig))
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

            if(mConfig.SpecificChrRegions != null && mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            if(mConfig.SpecificChrRegions != null && !mConfig.SpecificChrRegions.Regions.isEmpty())
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

    private void addCacheUnmappedReadFutures(final List<CompletableFuture<Void>> futures, final RemoteReadHandler remoteReadHandler,
            final ExecutorService executorService)
    {
        if(mConfig.SpecificChrRegions == null ||
                !mConfig.SpecificChrRegions.hasFilters() ||
                mConfig.SpecificChrRegions.Chromosomes.contains(CHR_UNMAPPED))
        {
            // write all unmapped reads to hash bams
            int numTasks = Math.max(mConfig.Threads / 10, 1);
            for(int i = 0; i < numTasks; ++i)
            {
                final int taskId = i;
                Runnable task = () -> remoteReadHandler.cacheAllUnmappedReads(numTasks, taskId);
                futures.add(CompletableFuture.runAsync(task, executorService));
            }
        }
    }

    // since each thread write to their own set of FASTQ files, we might need to merge them together
    private void mergeThreadFastqFiles(Collection<FastqWriterCache> fastqWriterCaches, ExecutorService executorService)
            throws ExecutionException, InterruptedException
    {
        // gather all the tasks
        final List<Runnable> tasks = new ArrayList<>();

        switch(mConfig.SplitMode)
        {
            case NONE:
            {
                BT_LOGGER.info("start merging fastqs");
                String filePrefix = mConfig.formFilePrefix("", "", true);
                String r1Fastq = formFilename(filePrefix, R1);
                String r2Fastq = formFilename(filePrefix, R2);
                String unpairedFastq = formFilename(filePrefix, UNPAIRED);

                // gather all the fastq writers
                final List<FastqWriter> fastqWriters = fastqWriterCaches.stream().map(FastqWriterCache::getFastqWriter)
                        .collect(Collectors.toList());

                tasks.add(() -> mergeFastqs(fastqWriters.stream().map(FastqWriter::getFastqR1), r1Fastq, false));
                tasks.add(() -> mergeFastqs(fastqWriters.stream().map(FastqWriter::getFastqR2), r2Fastq, false));
                tasks.add(() -> mergeFastqs(fastqWriters.stream().map(FastqWriter::getFastqUnpaired), unpairedFastq, true));
                break;
            }
            case READ_GROUP:
            {
                BT_LOGGER.info("start merging fastqs by read groups");
                for(SAMReadGroupRecord readGroup : ToFastqUtils.getReadGroups(mConfig))
                {
                    // we need to store all the read group details in the file name, but use slash-t instead of tab
                    // this is what BWA -R flag expects
                    String filePrefix = mConfig.OutputDir + readGroup.getSAMString().replace("\t", "\\t");
                    String r1Fastq = formFilename(filePrefix, R1);
                    String r2Fastq = formFilename(filePrefix, R2);
                    String unpairedFastq = formFilename(filePrefix, UNPAIRED);

                    BT_LOGGER.info("RG({}) R1 fastq({}) R2 fastq({})", readGroup.getId(), r1Fastq, r2Fastq);

                    // gather all the fastq writers
                    List<FastqWriter> fastqWriters = fastqWriterCaches.stream().map(o -> o.getReadGroupFastqWriter(readGroup.getId()))
                            .collect(Collectors.toList());

                    tasks.add(() -> mergeFastqs(fastqWriters.stream().map(FastqWriter::getFastqR1), r1Fastq, false));
                    tasks.add(() -> mergeFastqs(fastqWriters.stream().map(FastqWriter::getFastqR2), r2Fastq, false));
                    tasks.add(() -> mergeFastqs(fastqWriters.stream().map(FastqWriter::getFastqUnpaired), unpairedFastq, true));
                }
                break;
            }
            case THREAD:
                // do not need to do anything
        }

        if(!tasks.isEmpty())
        {
            final List<CompletableFuture<Void>> futures = new ArrayList<>();
            // submit all to thread pool
            tasks.forEach(o -> futures.add(CompletableFuture.runAsync(o, executorService)));
            // wait for completion
            CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();
            BT_LOGGER.info("finished merging fastqs");
        }
    }

    public static void main(final String[] args) throws ExecutionException, InterruptedException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        // set all thread exception handler
        // we must do this otherwise unhandled exception in other threads might not be reported
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            BT_LOGGER.fatal("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        BamToFastq bamToFastq = new BamToFastq(configBuilder);
        bamToFastq.run();
    }
}
