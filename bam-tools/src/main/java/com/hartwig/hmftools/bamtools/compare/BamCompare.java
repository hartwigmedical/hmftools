package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.executeRunnables;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamCompare
{
    private final CompareConfig mConfig;

    public BamCompare(final ConfigBuilder configBuilder)
    {
        mConfig = new CompareConfig(configBuilder);
    }

    public void run()
    {
        BT_LOGGER.info("starting bam comparison, writing output to {}", mConfig.OutputFile);

        long startTimeMs = System.currentTimeMillis();

        try(ReadWriter readWriter = new ReadWriter(mConfig))
        {
            if(!readWriter.initialised())
                System.exit(1);

            UnmatchedReadHandler unmatchedReadHandler = new UnmatchedReadHandler(mConfig);
            SamReaderFactory samReaderFactory = CompareUtils.makeSamReaderFactory(mConfig);
            List<Runnable> tasks = new ArrayList<>();
            List<PartitionReader> partitionTasks = new ArrayList<>();

            try(BamReaderProvider origThreadBamReader = BamReaderProvider.makeThreadLocal(samReaderFactory, new File(mConfig.OrigBamFile));
                    BamReaderProvider newThreadBamBReader = BamReaderProvider.makeThreadLocal(samReaderFactory, new File(mConfig.NewBamFile)))
            {
                // process unmapped
                addProcessUnmappedTasks(tasks, origThreadBamReader, newThreadBamBReader, unmatchedReadHandler);

                List<BamPartition> partitions = createPartitions();
                BT_LOGGER.info("splitting {} partitions across {} threads", partitions.size(), mConfig.Threads);

                for(BamPartition bamPartition : partitions)
                {
                    PartitionReader partitionReader = new PartitionReader(bamPartition.toString(), mConfig, bamPartition,
                            origThreadBamReader, newThreadBamBReader, readWriter, unmatchedReadHandler);
                    partitionTasks.add(partitionReader);
                    tasks.add(partitionReader);
                }

                if(!executeRunnables(tasks, mConfig.Threads))
                    System.exit(1);
            }

            Statistics combinedStats = new Statistics();
            partitionTasks.forEach(x -> combinedStats.merge(x.stats()));
            partitionTasks.clear();

            unmatchedReadHandler.closeHashBamWriters();

            List<BamReaderProvider> hashBamReaderProviders = new ArrayList<>();

            // now process all hash bam pairs
            for(Map.Entry<Integer, Pair<File, File>> hashBamPairEntries : unmatchedReadHandler.getHashBamPairs().entrySet())
            {
                Pair<File, File> hashBamPair = hashBamPairEntries.getValue();
                BamReaderProvider origHashBamReaderProvider = BamReaderProvider.of(samReaderFactory.open(hashBamPair.getLeft()));
                BamReaderProvider newHashBamReaderProvider = BamReaderProvider.of(samReaderFactory.open(hashBamPair.getRight()));
                partitionTasks.add(new PartitionReader(String.format("hash bam %d", hashBamPairEntries.getKey()),
                        mConfig, BamPartition.ofWholeBam(), origHashBamReaderProvider,
                        newHashBamReaderProvider, readWriter, null));
                hashBamReaderProviders.add(origHashBamReaderProvider);
                hashBamReaderProviders.add(newHashBamReaderProvider);
            }

            if(!executeRunnables(partitionTasks, mConfig.Threads))
                System.exit(1);

            hashBamReaderProviders.forEach(BamReaderProvider::close);

            partitionTasks.forEach(x -> combinedStats.merge(x.stats()));

            BT_LOGGER.printf(Level.INFO, "summary: reads(orig=%,d new=%,d) diffs(%,d)",
                    combinedStats.OrigReadCount, combinedStats.NewReadCount, combinedStats.DiffCount);
        }

        BT_LOGGER.info("BamCompare complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<BamPartition> createPartitions()
    {
        // get sequences of both bam files
        final List<SAMSequenceRecord> sequenceRecords = new ArrayList<>();
        SamReaderFactory samReaderFactory = CompareUtils.makeSamReaderFactory(mConfig);
        try(SamReader samReaderOrig = samReaderFactory.open(new File(mConfig.OrigBamFile));
            SamReader samReaderNew = samReaderFactory.open(new File(mConfig.NewBamFile)))
        {
            sequenceRecords.addAll(samReaderOrig.getFileHeader().getSequenceDictionary().getSequences());
            sequenceRecords.addAll(samReaderNew.getFileHeader().getSequenceDictionary().getSequences());
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        List<BamPartition> partitions = new ArrayList<>();

        Set<String> processedSequences = new HashSet<>();
        for(final SAMSequenceRecord sequenceRecord : sequenceRecords)
        {
            String chromosome = sequenceRecord.getSequenceName();

            if(!processedSequences.add(chromosome))
                continue;

            if(mConfig.SpecificChrRegions != null && mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            if(mConfig.SpecificChrRegions != null && !mConfig.SpecificChrRegions.Regions.isEmpty())
            {
                partitionChromosome(chromosome, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize)
                        .forEach(region -> partitions.add(BamPartition.ofRegion(region)));
            }
            else
            {
                buildPartitions(chromosome, sequenceRecord.getSequenceLength(), mConfig.PartitionSize)
                        .forEach(region -> partitions.add(BamPartition.ofRegion(region)));
            }
        }

        return partitions;
    }

    // pass all unmapped reads to the unmatched read handler
    private void addProcessUnmappedTasks(List<Runnable> tasks, BamReaderProvider origBamReaderProvider,
            BamReaderProvider newBamReaderProvider, UnmatchedReadHandler unmatchedReadHandler)
    {
        if(!mConfig.ignoreUnmapped())
        {
            tasks.add(() ->
            {
                try(SAMRecordIterator itr = origBamReaderProvider.getBamReader().queryUnmapped())
                {
                    long numReads = unmatchedReadHandler.handleOrigBamReads(itr);
                    BT_LOGGER.printf(Level.DEBUG, "finished writing %,d unmapped orig bam reads to hash bams", numReads);
                }
            });

            tasks.add(() ->
            {
                try(SAMRecordIterator itr = newBamReaderProvider.getBamReader().queryUnmapped())
                {
                    long numReads = unmatchedReadHandler.handleNewBamReads(itr);
                    BT_LOGGER.printf(Level.DEBUG, "finished writing %,d unmapped new bam reads to hash bams", numReads);
                }
            });
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        // set all thread exception handler
        // we must do this otherwise unhandled exception in other threads might not be reported
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            BT_LOGGER.fatal("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        BamCompare bamCompare = new BamCompare(configBuilder);
        bamCompare.run();
    }
}
