package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.bamtools.checker.PartitionThread.SORTED_BAM_ID;
import static com.hartwig.hmftools.bamtools.checker.PartitionThread.UNSORTED_BAM_ID;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.PartitionTask.splitRegionsIntoPartitions;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Queue;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.TaskQueue;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;

public class BamChecker
{
    private final CheckConfig mConfig;

    public BamChecker(final ConfigBuilder configBuilder)
    {
        mConfig = new CheckConfig(configBuilder);
    }

    public void run()
    {
        BT_LOGGER.info("starting BamChecker for file({})", mConfig.BamFile);

        long startTimeMs = System.currentTimeMillis();

        SAMFileHeader fileHeader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile)).getFileHeader();

        FragmentCache fragmentCache = new FragmentCache(fileHeader);

        List<PartitionThread> partitionThreads = createPartitionThreads(fragmentCache);
        List<Thread> allThreads = Lists.newArrayList(partitionThreads);

        if(!runThreadTasks(allThreads))
            System.exit(1);

        BT_LOGGER.info("all partition tasks complete, mins({})", runTimeMinsStr(startTimeMs));

        fragmentCache.logFinalStats();

        FragmentStats combinedStats = new FragmentStats();
        combinedStats.merge(fragmentCache.stats());
        partitionThreads.forEach(x -> combinedStats.merge(x.stats()));

        BT_LOGGER.debug("total stats: {}", combinedStats.toString());

        List<SAMRecord> incompleteReads = fragmentCache.extractReads();

        if(!incompleteReads.isEmpty())
        {
            BT_LOGGER.debug("incomplete {} reads", incompleteReads.size());

            long primaryCount = incompleteReads.stream().filter(x -> !x.getSupplementaryAlignmentFlag()).count();

            if(primaryCount > 0)
            {
                BT_LOGGER.warn("incomplete {} primary reads", primaryCount);
            }

            if(mConfig.WriteIncompleteFragments)
            {
                writeIncompleteReads(incompleteReads);
            }

            fragmentCache.clear();
        }

        if(mConfig.writeBam())
        {
            finaliseBam(partitionThreads, incompleteReads);
        }

        BT_LOGGER.info("BamChecker complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<PartitionThread> createPartitionThreads(final FragmentCache fragmentCache)
    {
        List<PartitionThread> partitionThreads = Lists.newArrayListWithCapacity(mConfig.Threads);

        List<ChrBaseRegion> partitionRegions = PartitionThread.splitRegionsIntoPartitions(mConfig);

        if(partitionRegions.isEmpty())
            return Collections.emptyList();

        BT_LOGGER.debug("splitting {} partition regions across {} threads", partitionRegions.size(), mConfig.Threads);

        Queue<ChrBaseRegion> partitionQueue = new ConcurrentLinkedQueue<>();
        partitionQueue.addAll(partitionRegions);

        TaskQueue taskQueue = new TaskQueue(partitionQueue, "partitions", 0); // log on completed partitions

        for(int i = 0; i < mConfig.Threads; ++i)
        {
            PartitionThread partitionThread = new PartitionThread(mConfig, fragmentCache, taskQueue, i);
            partitionThreads.add(partitionThread);
        }

        return partitionThreads;
    }

    private void finaliseBam(final List<PartitionThread> partitionThreads, final List<SAMRecord> incompleteReads)
    {
        if(!incompleteReads.isEmpty() && !mConfig.DropIncompleteFragments)
        {
            partitionThreads.get(0).writeIncompleteReads(incompleteReads);
        }

        // write fully unmapped reads
        if(!mConfig.SkipUnmapped)
        {
            BT_LOGGER.debug("writing unmapped reads");

            PartitionThread writerThread = partitionThreads.size() > 1 ? partitionThreads.get(1) : partitionThreads.get(0);
            writerThread.writeUnmappedReads();
        }

        partitionThreads.forEach(x -> x.close());

        // sort and merge the interim BAMs
        BT_LOGGER.debug("sorting {} thread BAMs", partitionThreads.size());

        BamToolName toolName = fromPath(mConfig.BamToolPath);

        List<String> unsortedBams = Lists.newArrayListWithCapacity(mConfig.Threads);
        List<String> sortedBams = Lists.newArrayListWithCapacity(mConfig.Threads);
        List<BamSortTask> sortTasks = Lists.newArrayListWithCapacity(mConfig.Threads);

        for(PartitionThread partitionThread : partitionThreads)
        {
            String unsortedBamFilename = partitionThread.bamFilename();
            unsortedBams.add(unsortedBamFilename);

            String sortedBamFilename = unsortedBamFilename.replaceAll(UNSORTED_BAM_ID, SORTED_BAM_ID);
            sortedBams.add(sortedBamFilename);

            sortTasks.add(new BamSortTask(unsortedBamFilename, sortedBamFilename));
        }

        List<Callable> threadTasks = sortTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(threadTasks, mConfig.Threads))
        {
            BT_LOGGER.error("failed to sort {} BAMs", sortedBams.size());
            System.exit(1);
        }

        // merge sorted BAMs
        String finalBam = mConfig.formFilename("final", BAM_EXTENSION);
        if(!BamOperations.mergeBams(toolName, mConfig.BamToolPath, finalBam, sortedBams, mConfig.Threads))
        {
            BT_LOGGER.error("error merging sorted BAMs");
            System.exit(1);
        }

        // index final BAM
        String finalBamIndex = finalBam + BAM_INDEX_EXTENSION;
        if(!Files.exists(Paths.get(finalBamIndex)))
        {
            BamOperations.indexBam(toolName, mConfig.BamToolPath, finalBam, mConfig.Threads);
        }

        // clean-up interim BAMs
        try
        {
            for(int i = 0; i < sortedBams.size(); ++i)
            {
                Files.deleteIfExists(Paths.get(unsortedBams.get(i)));

                String sortedBamFilename = sortedBams.get(i);
                Files.deleteIfExists(Paths.get(sortedBamFilename));
                Files.deleteIfExists(Paths.get(sortedBamFilename + BAM_INDEX_EXTENSION));
            }
        }
        catch(IOException e)
        {
            BT_LOGGER.error("error deleting interim file: {}", e.toString());
            System.exit(1);
        }
    }

    private class BamSortTask implements Callable
    {
        private final String mInputBam;
        private final String mOutputBam;

        public BamSortTask(final String inputBam, final String outputBam)
        {
            mInputBam = inputBam;
            mOutputBam = outputBam;
        }

        @Override
        public Long call()
        {
            BamToolName toolName = fromPath(mConfig.BamToolPath);
            BamOperations.sortBam(toolName, mConfig.BamToolPath, mInputBam, mOutputBam, 1);
            return (long)0;
        }
    }

    private void writeIncompleteReads(final List<SAMRecord> incompleteReads)
    {
        try
        {
            String filename = mConfig.formFilename("incomplete_read", TSV_EXTENSION);
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("ReadId").add("Chromosome").add("PosStart").add("PosEnd").add("Cigar");
            sj.add("MateChr").add("MatePosStart").add("SuppData").add("Flags");
            sj.add("FirstInPair").add("Unmapped").add("MateUnmapped").add("Supplementary");

            writer.write(sj.toString());

            writer.newLine();

            for(SAMRecord read : incompleteReads)
            {
                sj = new StringJoiner(TSV_DELIM);

                sj.add(read.getReadName());
                sj.add(read.getContig());
                sj.add(String.valueOf(read.getAlignmentStart()));
                sj.add(String.valueOf(read.getAlignmentEnd()));
                sj.add(read.getCigarString());

                List<SupplementaryReadData> suppData = SupplementaryReadData.extractAlignments(read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

                sj.add(read.getMateReferenceName());
                sj.add(String.valueOf(read.getMateAlignmentStart()));

                if(suppData != null)
                {
                    String suppDataStr = suppData.stream().map(x -> x.asDelimStr()).collect(Collectors.joining("|"));
                    sj.add(suppDataStr);
                }
                else
                {
                    sj.add("N/A");
                }

                sj.add(String.valueOf(read.getFlags()));

                sj.add(String.valueOf(read.getFirstOfPairFlag()));
                sj.add(String.valueOf(read.getReadUnmappedFlag()));
                sj.add(String.valueOf(read.getMateUnmappedFlag()));
                sj.add(String.valueOf(read.getSupplementaryAlignmentFlag()));

                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to writer incomplete read writer: {}", e.toString());
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CheckConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BamChecker bamChecker = new BamChecker(configBuilder);
        bamChecker.run();
    }
}
