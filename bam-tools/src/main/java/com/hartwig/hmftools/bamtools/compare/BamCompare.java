package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.PartitionReader.NAME_UNMAPPED;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
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
        BT_LOGGER.info("starting BAM comparison, writing output to {}", mConfig.OutputFile);

        long startTimeMs = System.currentTimeMillis();

        ReadWriter readWriter = new ReadWriter(mConfig);

        if(!readWriter.initialised())
            System.exit(1);

        UnmatchedReadHandler unmatchedReadHandler = new UnmatchedReadHandler(mConfig);

        List<PartitionReader> partitionTasks = new ArrayList<>();

        File origBamFile = new File(mConfig.OrigBamFile);
        File newBamFile = new File(mConfig.NewBamFile);

        List<ChrBaseRegion> partitions = createPartitions();
        BT_LOGGER.info("splitting {} partitions across {} threads", partitions.size(), mConfig.Threads);

        for(ChrBaseRegion partition : partitions)
        {
            PartitionReader partitionReader = new PartitionReader(
                    partition.toString(), mConfig, partition, origBamFile, newBamFile, readWriter, unmatchedReadHandler);

            partitionTasks.add(partitionReader);
        }

        if(!mConfig.ignoreUnmapped())
        {
            PartitionReader partitionReader = new PartitionReader(
                    NAME_UNMAPPED, mConfig, null, origBamFile, newBamFile, readWriter, unmatchedReadHandler);

            partitionTasks.add(partitionReader);
        }

        BT_LOGGER.info("partition processing complete");

        List<Callable<Void>> callableList = partitionTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        Statistics combinedStats = new Statistics();
        partitionTasks.forEach(x -> combinedStats.merge(x.stats()));
        partitionTasks.clear();

        unmatchedReadHandler.closeHashBamWriters();

        if(!unmatchedReadHandler.getHashBamPairs().isEmpty())
        {
            BT_LOGGER.info("processing cached hash BAMs");

            // now process all hash bam pairs
            for(Map.Entry<Integer, Pair<File, File>> hashBamPairEntries : unmatchedReadHandler.getHashBamPairs().entrySet())
            {
                Pair<File,File> hashBamPair = hashBamPairEntries.getValue();

                partitionTasks.add(new PartitionReader(String.format("hash bam %d", hashBamPairEntries.getKey()),
                        mConfig, null, hashBamPair.getLeft(), hashBamPair.getRight(), readWriter, null));
           }

            callableList = partitionTasks.stream().collect(Collectors.toList());

            if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
                System.exit(1);

            partitionTasks.forEach(x -> combinedStats.merge(x.stats()));
        }

        readWriter.close();

        BT_LOGGER.printf(Level.INFO, "summary: reads(orig=%,d new=%,d) diffs(%,d)",
                combinedStats.OrigReadCount, combinedStats.NewReadCount, combinedStats.DiffCount);

        BT_LOGGER.info("BamCompare complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<ChrBaseRegion> createPartitions()
    {
        // get sequences of both BAM headers
        final List<SAMSequenceRecord> sequenceRecords = new ArrayList<>();
        SamReaderFactory samReaderFactory = CompareUtils.makeSamReaderFactory(mConfig);

        SamReader samReaderOrig = samReaderFactory.open(new File(mConfig.OrigBamFile));
        List<SAMSequenceRecord> sequencesOriginal = samReaderOrig.getFileHeader().getSequenceDictionary().getSequences();
        sequenceRecords.addAll(sequencesOriginal);

        SamReader samReaderNew = samReaderFactory.open(new File(mConfig.NewBamFile));
        List<SAMSequenceRecord> sequencesNew = samReaderNew.getFileHeader().getSequenceDictionary().getSequences();
        sequenceRecords.addAll(sequencesNew);

        List<ChrBaseRegion> partitions = Lists.newArrayList();

        Set<String> processedSequences = new HashSet<>();

        for(SAMSequenceRecord sequenceRecord : sequenceRecords)
        {
            String chromosome = sequenceRecord.getSequenceName();

            if(!processedSequences.add(chromosome))
                continue;

            if(mConfig.SpecificChrRegions != null && mConfig.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            if(mConfig.SpecificChrRegions != null && !mConfig.SpecificChrRegions.Regions.isEmpty())
            {
                partitions.addAll(partitionChromosome(
                        chromosome, mConfig.RefGenVersion, mConfig.SpecificChrRegions.Regions, mConfig.PartitionSize));
            }
            else
            {
                partitions.addAll(buildPartitions(chromosome, sequenceRecord.getSequenceLength(), mConfig.PartitionSize));
            }
        }

        return partitions;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CompareConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BamCompare bamCompare = new BamCompare(configBuilder);
        bamCompare.run();
    }
}
