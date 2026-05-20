package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.CompareConfig.FULLY_UNMAPPED_PARTITION;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.NotNull;

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

        final IgnoreRegionIndex ignoreRegionIndex = mConfig.IgnoreRegionFiles.isEmpty()
                ? null : IgnoreRegionIndex.load(mConfig.IgnoreRegionFiles);
        final Set<String> ignoredReadNames = ignoreRegionIndex != null
                ? Collections.newSetFromMap(new ConcurrentHashMap<>()) : null;

        List<PartitionReader> partitionTasks = new ArrayList<>();

        File origBamFile = new File(mConfig.OrigBamFile);
        File newBamFile = new File(mConfig.NewBamFile);

        List<ChrBaseRegion> partitions = createPartitions();

        for(ChrBaseRegion partition : partitions)
        {
            PartitionReader partitionReader = new PartitionReader(
                    partition.toString(), mConfig, partition, origBamFile, newBamFile, readWriter, unmatchedReadHandler,
                    ignoreRegionIndex, ignoredReadNames);

            partitionTasks.add(partitionReader);
        }

        boolean includeFullyUnmapped = !mConfig.StandardChromosomes
                && (!mConfig.SpecificChrRegions.hasFilters() || mConfig.SpecificChrRegions.Chromosomes.contains(FULLY_UNMAPPED_PARTITION));

        if(includeFullyUnmapped)
        {
            PartitionReader partitionReader = new PartitionReader(
                    FULLY_UNMAPPED_PARTITION, mConfig, null, origBamFile, newBamFile, readWriter, unmatchedReadHandler,
                    ignoreRegionIndex, ignoredReadNames);

            partitionTasks.add(partitionReader);
        }

        BT_LOGGER.info("splitting {} partitions across {} threads", partitionTasks.size(), mConfig.Threads);

        // pre-pass: scan both BAMs over each partition and collect readnames whose any alignment overlaps
        // an ignore region by > 50% of its reference span. Populates the shared ignoredReadNames set used
        // by excludeRead in the compare wave below.
        if(ignoreRegionIndex != null)
        {
            long scanStartMs = System.currentTimeMillis();
            BT_LOGGER.info("ignore-region pre-pass: scanning {} partitions for readnames to drop", partitionTasks.size());

            List<Callable<Void>> scanCallables = partitionTasks.stream()
                    .map(p -> (Callable<Void>) () -> { p.scanIgnoredReadNames(); return null; })
                    .collect(Collectors.toList());

            if(!TaskExecutor.executeTasks(scanCallables, mConfig.Threads))
                System.exit(1);

            long scanSecs = (System.currentTimeMillis() - scanStartMs) / 1000;
            BT_LOGGER.info("ignore-region pre-pass complete: {} readname(s) flagged, {} sec",
                    ignoredReadNames.size(), scanSecs);
        }

        List<Callable<Void>> callableList = partitionTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        BT_LOGGER.info("partition processing complete");

        Statistics combinedStats = new Statistics();
        partitionTasks.forEach(x -> combinedStats.merge(x.stats()));
        partitionTasks.clear();

        unmatchedReadHandler.closeHashBamWriters();

        System.gc(); // probably unnecessary

        if(!unmatchedReadHandler.getHashBamPairs().isEmpty())
        {
            BT_LOGGER.info("processing cached hash BAMs");

            // now process all hash bam pairs
            for(Map.Entry<Integer, Pair<File, File>> hashBamPairEntries : unmatchedReadHandler.getHashBamPairs().entrySet())
            {
                Pair<File,File> hashBamPair = hashBamPairEntries.getValue();

                partitionTasks.add(new PartitionReader(String.format("hash bam %d", hashBamPairEntries.getKey()),
                        mConfig, null, hashBamPair.getLeft(), hashBamPair.getRight(), readWriter, null,
                        ignoreRegionIndex, ignoredReadNames));
           }

            callableList = partitionTasks.stream().collect(Collectors.toList());

            if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
                System.exit(1);

            partitionTasks.forEach(x -> combinedStats.merge(x.stats()));
        }

        readWriter.close();

        BT_LOGGER.printf(Level.INFO, "summary: reads(orig=%,d new=%,d) matched(%,d) diffs(%,d) ignored(%,d)",
                combinedStats.OrigReadCount, combinedStats.NewReadCount, combinedStats.Matched, combinedStats.DiffCount,
                combinedStats.IgnoredReadRecords);

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
