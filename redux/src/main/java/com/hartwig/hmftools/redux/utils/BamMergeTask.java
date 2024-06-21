package com.hartwig.hmftools.redux.utils;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.utils.BamMerger.formSequenceBamFilename;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamMergeTask extends Thread
{
    private final List<String> mInputBams;
    private final Queue<SAMSequenceRecord> mSAMSequences;
    private final String mOutputBamPrefix;

    private final List<BamSequenceReader> mActiveBamReaders;
    private final List<BamSequenceReader> mFinishedBamReaders;
    private final String mRefGenomeFile;
    private SAMFileWriter mSamFileWriter;

    private int mReorderCount;

    public BamMergeTask(
            final List<String> inputBams, final String refGenomeFile, final Queue<SAMSequenceRecord> sequences, final String outputBamPrefix)
    {
        mSAMSequences = sequences;
        mOutputBamPrefix = outputBamPrefix;
        mInputBams = inputBams;
        mRefGenomeFile = refGenomeFile;

        mActiveBamReaders = Lists.newArrayListWithCapacity(inputBams.size());
        mFinishedBamReaders = Lists.newArrayListWithCapacity(inputBams.size());

        mSamFileWriter = null;
        mReorderCount = 0;

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                SAMSequenceRecord sequence = mSAMSequences.remove();
                processSequence(sequence);
            }
            catch(NoSuchElementException e)
            {
                RD_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private static final int LOG_COUNT = 10_000_000;

    private void processSequence(final SAMSequenceRecord sequence)
    {
        // open and prepare each BAM, adding in start order
        for(String bamFile : mInputBams)
        {
            BamSequenceReader bamReader = new BamSequenceReader(mRefGenomeFile, bamFile, sequence);

            if(bamReader.finished())
            {
                mFinishedBamReaders.add(bamReader);
            }
            else
            {
                addBamReaderInPosition(bamReader);
            }
        }

        if(mActiveBamReaders.isEmpty())
        {
            RD_LOGGER.debug("contig({}) no BAM files with records found", sequence.getSequenceName(), mActiveBamReaders.size());
            return;
        }

        String outputBam = formSequenceBamFilename(mOutputBamPrefix, sequence);
        mSamFileWriter = initialiseWriter(outputBam);

        RD_LOGGER.debug("contig({}) merging {} BAMs", sequence.getSequenceName(), mActiveBamReaders.size());

        // begin the merge process:
        // 1. take the top record from the first BAM reader
        // 2. check the next record vs the next top record(s), and if not still top then move the BAM reader to its new position
        // 3. repeat the process
        BamSequenceReader topWriter = mActiveBamReaders.get(0);

        long recordCount = 0;

        while(!mActiveBamReaders.isEmpty())
        {
            mSamFileWriter.addAlignment(topWriter.current());
            ++recordCount;

            if((recordCount % LOG_COUNT) == 0)
            {
                RD_LOGGER.trace("contig({}) merged {} records, readers(active={} finished={}) reorders({})",
                        sequence.getSequenceName(), recordCount, mActiveBamReaders.size(), mFinishedBamReaders.size(), mReorderCount);
            }

            topWriter.moveNext();

            // test most likely scnenario first
            if(topWriter.current() != null && mActiveBamReaders.size() > 1 && topWriter.isLowerOrEqualWith(mActiveBamReaders.get(1)))
            {
                continue;
            }

            mActiveBamReaders.remove(0);

            if(topWriter.finished())
            {
                RD_LOGGER.trace("contig({}) bam({}) finished", sequence.getSequenceName(), topWriter.filename());
                mFinishedBamReaders.add(topWriter);
            }
            else
            {
                ++mReorderCount;

                // move the top reader to its new position in the list and repeat the process
                addBamReaderInPosition(topWriter);
            }

            if(mActiveBamReaders.isEmpty())
                break;

            topWriter = mActiveBamReaders.get(0);
        }

        RD_LOGGER.debug("contig({}) merged {} BAM files with {} records, reorder count({})",
                sequence.getSequenceName(), mInputBams.size(), recordCount, mReorderCount);

        mSamFileWriter.close();
    }

    private void addBamReaderInPosition(final BamSequenceReader bamReader)
    {
        int index = 0;
        while(index < mActiveBamReaders.size())
        {
            if(bamReader.isHigherThan(mActiveBamReaders.get(index)))
                ++index;
            else
                break;
        }

        mActiveBamReaders.add(index, bamReader);
    }

    private SAMFileWriter initialiseWriter(final String outputBam)
    {
        String sampleBam = mInputBams.get(0);
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mRefGenomeFile)).open(new File(sampleBam));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        // need to check this - must be unsorted to write as
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, true, new File(outputBam));
    }
}
