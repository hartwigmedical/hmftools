package com.hartwig.hmftools.common.bamops;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bamops.BamOperations.BOP_LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CRAM_INDEX_EXTENSION;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamMerger
{
    private final String mOutputBamPrefix;
    private final List<String> mInputBams;
    private final String mRefGenomeFile;
    private final String mBamToolPath;
    private final int mThreads;
    private final boolean mKeepInterimBams;

    protected static final String UNMAPPED_READS = "unmapped";

    public BamMerger(
            final String outputBam, final List<String> inputBams, final String refGenomeFile, final String bamToolPath,
            final int threads, final boolean keepInterimBams)
    {
        String outputBamPrefix = outputBam.substring(0, outputBam.indexOf(BAM_EXTENSION));
        mOutputBamPrefix = outputBamPrefix;
        mBamToolPath = bamToolPath;
        mInputBams = inputBams;
        mRefGenomeFile = refGenomeFile;
        mThreads = max(threads, 1);
        mKeepInterimBams = keepInterimBams;
    }

    public int inputBamCount() { return mInputBams.size(); }

    public boolean merge()
    {
        if(mInputBams.size() < 2)
        {
            BOP_LOGGER.warn("invalid input BAM file count({})", mInputBams.size());
            return false;
        }

        buildIndexFiles();

        // collect all unique sequences from the BAMs
        RefGenomeSource refGenome = loadRefGenome(mRefGenomeFile);
        List<SAMSequenceRecord> sequences = refGenome.refGenomeFile().getSequenceDictionary().getSequences();

        int sequenceInfoCount = mThreads;
        List<SequenceInfo> sequenceIntervals = formSequenceIntervals(sequences, mOutputBamPrefix, sequenceInfoCount);

        Queue<SequenceInfo> sequenceIntervalsQueue = new ConcurrentLinkedQueue<>();

        sequenceIntervals.forEach(x -> sequenceIntervalsQueue.add(x));

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(sequences.size(), mThreads); ++i)
        {
            BamMergeTask bamMergeTask = new BamMergeTask(mInputBams, mRefGenomeFile, sequenceIntervalsQueue);
            workers.add(bamMergeTask);
        }

        UnmappedMergeTask unmappedMergeTask = new UnmappedMergeTask(mInputBams, mRefGenomeFile, mOutputBamPrefix);
        workers.add(unmappedMergeTask);

        BOP_LOGGER.debug("splitting {} sequence merges across {} threads", sequences.size(), mThreads);

        if(!runThreadTasks(workers))
            System.exit(1);

        BOP_LOGGER.debug("all sequence merge tasks complete");

        makeFinalBam(sequenceIntervals);

        BOP_LOGGER.debug("BAM merge complete");

        return true;
    }

    public static List<SequenceInfo> formSequenceIntervals(
            final List<SAMSequenceRecord> sequences, final String bamPrefix, final int sequenceInfoCount)
    {
        long totalLength = sequences.stream().mapToLong(x -> x.getSequenceLength()).sum();
        long intervalLength = (int)ceil(totalLength / (double)sequenceInfoCount);

        List<SequenceInfo> sequenceInfoItems = Lists.newArrayListWithCapacity(sequenceInfoCount);

        long currentLength = 0;
        int nextSequenceStart = 1;
        SequenceInfo currentSeqInfo = null;

        for(SAMSequenceRecord sequence : sequences)
        {
            nextSequenceStart = 1;
            int remainingSequenceLength = sequence.getSequenceLength() - nextSequenceStart + 1;

            if(currentSeqInfo == null)
            {
                String seqBam = formBamFilename(bamPrefix, String.valueOf(sequenceInfoItems.size()));
                currentSeqInfo = new SequenceInfo(sequenceInfoItems.size(), seqBam);
            }

            while(currentLength + remainingSequenceLength >= intervalLength)
            {
                int remainingLength = (int) (intervalLength - currentLength);
                int sequenceEnd = nextSequenceStart + remainingLength - 1;

                currentSeqInfo.Intervals.add(new QueryInterval(sequence.getSequenceIndex(), nextSequenceStart, sequenceEnd));

                if(!sequenceInfoItems.contains(currentSeqInfo))
                    sequenceInfoItems.add(currentSeqInfo);

                // start a new sequence info
                String seqBam = formBamFilename(bamPrefix, String.valueOf(sequenceInfoItems.size()));
                currentSeqInfo = new SequenceInfo(sequenceInfoItems.size(), seqBam);
                currentLength = 0;

                nextSequenceStart = sequenceEnd + 1;
                remainingSequenceLength = sequence.getSequenceLength() - nextSequenceStart + 1;
            }

            if(remainingSequenceLength <= 0)
                continue;

            currentLength += remainingSequenceLength;

            currentSeqInfo.Intervals.add(new QueryInterval(
                    sequence.getSequenceIndex(), nextSequenceStart, sequence.getSequenceLength()));

            if(!sequenceInfoItems.contains(currentSeqInfo))
                sequenceInfoItems.add(currentSeqInfo);
        }

        return sequenceInfoItems;
    }

    private void buildIndexFiles()
    {
        if(mBamToolPath == null)
            return;

        List<String> bamMissingIndexFiles = Lists.newArrayList();

        for(String inputBam : mInputBams)
        {
            String fileExtension = inputBam.endsWith(BAM_EXTENSION) ? BAM_INDEX_EXTENSION : CRAM_INDEX_EXTENSION;
            String indexFile = inputBam + fileExtension;

            if(!Files.exists(Paths.get(indexFile)))
                bamMissingIndexFiles.add(inputBam);
        }

        if(bamMissingIndexFiles.isEmpty())
            return;

        BOP_LOGGER.debug("building index files for {} files", bamMissingIndexFiles.size());

        BamToolName bamToolName = BamToolName.fromPath(mBamToolPath);

        for(String inputBam : bamMissingIndexFiles)
        {
            if(!BamOperations.indexBam(bamToolName, mBamToolPath, inputBam, mThreads))
                System.exit(1);
        }
    }

    protected static String formBamFilename(final String outputBamPrefix, final String bamFileId)
    {
        return outputBamPrefix + "_seq" + bamFileId + BAM_EXTENSION;
    }

    private void makeFinalBam(final List<SequenceInfo> sequenceIntervals)
    {
        String finalBam = mOutputBamPrefix + BAM_EXTENSION;

        List<String> interimBams = Lists.newArrayList();

        for(SequenceInfo sequenceInfo : sequenceIntervals)
        {
            if(Files.exists(Paths.get(sequenceInfo.BamFile)))
            {
                interimBams.add(sequenceInfo.BamFile);
            }
        }

        // add BAM with unmapped reads
        String unmappedBam = formBamFilename(mOutputBamPrefix, UNMAPPED_READS);
        interimBams.add(unmappedBam);

        BamToolName bamToolName = BamToolName.fromPath(mBamToolPath);
        BamOperations.concatenateBams(bamToolName, mBamToolPath, finalBam, interimBams, mThreads);

        if(!BamOperations.indexBam(bamToolName, mBamToolPath, finalBam, mThreads))
            System.exit(1);

        if(!mKeepInterimBams)
        {
            // clean-up interim files
            for(String interimBam : interimBams)
            {
                try
                {
                    Files.deleteIfExists(Paths.get(interimBam));
                    Files.deleteIfExists(Paths.get(interimBam + BAM_INDEX_EXTENSION));
                }
                catch(Exception e) {}
            }
        }
    }

    public static SAMFileHeader buildCombinedHeader(final List<String> bamFiles, final String refGenomeFile)
    {
        if(bamFiles.isEmpty())
            return null;

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile))
                .open(new File(bamFiles.get(0)));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        if(bamFiles.size() > 1)
        {
            for(int i = 1; i < bamFiles.size(); ++i)
            {
                SamReader nextReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile))
                        .open(new File(bamFiles.get(i)));

                for(SAMReadGroupRecord readGroupRecord : nextReader.getFileHeader().getReadGroups())
                {
                    if(fileHeader.getReadGroups().stream().noneMatch(x -> x.getReadGroupId().equals(readGroupRecord.getReadGroupId())))
                        fileHeader.addReadGroup(readGroupRecord);
                }

                if(!nextReader.getFileHeader().getProgramRecords().isEmpty())
                {
                    final SAMProgramRecord nextProgramRecord = nextReader.getFileHeader().getProgramRecords().get(0);
                    String newProgramId = String.format("%s.%d", nextProgramRecord.getId(), i);

                    if(fileHeader.getProgramRecords().stream().noneMatch(x -> x.getProgramGroupId().equals(newProgramId)))
                        fileHeader.addProgramRecord(new SAMProgramRecord(newProgramId, nextProgramRecord));
                }
            }
        }

        return fileHeader;
    }
}
