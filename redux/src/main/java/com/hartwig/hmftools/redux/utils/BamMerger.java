package com.hartwig.hmftools.redux.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamOperations;
import com.hartwig.hmftools.common.bam.BamToolName;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordIterator;
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

    protected static final String UNMAPPED_READS = "unmapped";

    public BamMerger(
            final String outputBam, final List<String> inputBams, final String refGenomeFile, final String bamToolPath,
            final int threads)
    {
        String outputBamPrefix = outputBam.substring(0, outputBam.indexOf(BAM_EXTENSION));
        mOutputBamPrefix = outputBamPrefix;
        mBamToolPath = bamToolPath;
        mInputBams = inputBams;
        mRefGenomeFile = refGenomeFile;
        mThreads = threads;
    }

    public boolean merge()
    {
        if(mInputBams.isEmpty())
            return false;

        // check input BAMs exist
        boolean hasMissing = false;

        for(String inputBam : mInputBams)
        {
            if(!Files.exists(Paths.get(inputBam)))
            {
                RD_LOGGER.error("missing input BAM: {}", inputBam);
                hasMissing = true;
            }
        }

        if(hasMissing)
            System.exit(1);

        buildIndexFiles();

        // collect all unique sequences from the BAMs
        RefGenomeSource refGenome = loadRefGenome(mRefGenomeFile);
        List<SAMSequenceRecord> sequences = refGenome.refGenomeFile().getSequenceDictionary().getSequences();

        Queue<SAMSequenceRecord> sequenceQueue = new ConcurrentLinkedQueue<>();

        sequences.forEach(x -> sequenceQueue.add(x));

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < min(sequences.size(), mThreads); ++i)
        {
            BamMergeTask bamMergeTask = new BamMergeTask(mInputBams, mRefGenomeFile, sequenceQueue, mOutputBamPrefix);
            workers.add(bamMergeTask);
        }

        UnmappedMergeTask unmappedMergeTask = new UnmappedMergeTask(mInputBams, mRefGenomeFile, mOutputBamPrefix);
        workers.add(unmappedMergeTask);

        RD_LOGGER.debug("splitting {} sequence merges across {} threads", sequences.size(), mThreads);

        if(!runThreadTasks(workers))
            System.exit(1);

        RD_LOGGER.info("all sequence merge tasks complete");

        // now sequentially add the sequence BAMs together
        // writeFinalBam(sequences);
        makeFinalBam(sequences);

        RD_LOGGER.info("BAM merge complete");

        return true;
    }

    private void buildIndexFiles()
    {
        if(mBamToolPath == null)
            return;

        List<String> bamMissingIndexFiles = Lists.newArrayList();

        for(String inputBam : mInputBams)
        {
            String indexFile = inputBam + BAM_INDEX_EXTENSION;

            if(!Files.exists(Paths.get(indexFile)))
                bamMissingIndexFiles.add(inputBam);
        }

        if(bamMissingIndexFiles.isEmpty())
            return;

        RD_LOGGER.info("building index files for {} files", bamMissingIndexFiles.size());

        BamToolName bamToolName = BamToolName.fromPath(mBamToolPath);

        for(String inputBam : bamMissingIndexFiles)
        {
            if(!BamOperations.indexBam(bamToolName, mBamToolPath, inputBam, mThreads))
                System.exit(1);
        }
    }

    protected static String formSequenceBamFilename(final String outputBamPrefix, final SAMSequenceRecord sequence)
    {
        return formBamFilename(outputBamPrefix, sequence.getSequenceName());
    }

    protected static String formBamFilename(final String outputBamPrefix, final String bamFileId)
    {
        return outputBamPrefix + "_seq" + bamFileId + BAM_EXTENSION;
    }

    private void makeFinalBam(final List<SAMSequenceRecord> sequences)
    {
        String finalBam = mOutputBamPrefix + BAM_EXTENSION;

        List<String> interimBams = Lists.newArrayList();

        for(SAMSequenceRecord sequence : sequences)
        {
            String sequenceBam = formSequenceBamFilename(mOutputBamPrefix, sequence);

            if(Files.exists(Paths.get(sequenceBam)))
            {
                interimBams.add(sequenceBam);
            }
        }

        // add BAM with unmapped reads
        String unmappedBam = formBamFilename(mOutputBamPrefix, UNMAPPED_READS);
        interimBams.add(unmappedBam);

        BamToolName bamToolName = BamToolName.fromPath(mBamToolPath);
        BamOperations.concatenateBams(bamToolName, mBamToolPath, finalBam, interimBams, mThreads);

        if(!BamOperations.indexBam(bamToolName, mBamToolPath, finalBam, mThreads))
            System.exit(1);

        // clean-up interim files
        for(String interimBam : interimBams)
        {
            try { Files.deleteIfExists(Paths.get(interimBam)); } catch(Exception e) {}
        }
    }

    private void writeFinalBam(final List<SAMSequenceRecord> sequences)
    {
        String finalBam = mOutputBamPrefix + BAM_EXTENSION;

        String sampleBam = mInputBams.get(0);
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mRefGenomeFile)).open(new File(sampleBam));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        // need to check this - must be unsorted to write as
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        SAMFileWriter bamWriter = new SAMFileWriterFactory().makeBAMWriter(fileHeader, true, new File(finalBam));

        for(SAMSequenceRecord sequence : sequences)
        {
            String sequenceBam = formSequenceBamFilename(mOutputBamPrefix, sequence);

            if(!Files.exists(Paths.get(sequenceBam)))
                continue;

            /*
            SamReader bamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mRefGenomeFile)).open(new File(sequenceBam));
            SAMRecordIterator iterator = bamReader.iterator();

            while(iterator.hasNext())
            {
                bamWriter.addAlignment(iterator.next());
            }
            */
        }

        RD_LOGGER.debug("writing unmapped reads");

        // finally add in unmapped reads
        for(String inputBam : mInputBams)
        {
            SamReader bamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mRefGenomeFile)).open(new File(inputBam));

            SAMRecordIterator iterator = bamReader.queryUnmapped();
            while(iterator.hasNext())
            {
                bamWriter.addAlignment(iterator.next());
            }
        }

        bamWriter.close();

        if(mBamToolPath != null)
        {
            BamToolName bamToolName = BamToolName.fromPath(mBamToolPath);
            BamOperations.indexBam(bamToolName, mBamToolPath, finalBam, mThreads);
        }
    }
}
