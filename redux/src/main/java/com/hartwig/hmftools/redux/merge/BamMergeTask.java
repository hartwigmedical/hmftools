package com.hartwig.hmftools.redux.merge;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamMergeTask extends Thread
{
    private final List<String> mInputBams;
    private final Queue<SequenceInfo> mSAMSequences;

    private final List<BamSequenceReader> mActiveBamReaders;
    private final List<BamSequenceReader> mFinishedBamReaders;
    private final String mRefGenomeFile;
    private SAMFileWriter mSamFileWriter;

    private int mReorderCount;

    public BamMergeTask(
            final List<String> inputBams, final String refGenomeFile, final Queue<SequenceInfo> sequences, final String outputBamPrefix)
    {
        mSAMSequences = sequences;
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
                SequenceInfo sequenceInfo = mSAMSequences.remove();
                processSequence(sequenceInfo);
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

    private void processSequence(final SequenceInfo sequenceInfo)
    {
        // open and prepare each BAM, adding in start order
        for(String bamFile : mInputBams)
        {
            BamSequenceReader bamReader = new BamSequenceReader(mRefGenomeFile, bamFile, sequenceInfo);

            if(bamReader.finished())
            {
                mFinishedBamReaders.add(bamReader);
            }
            else
            {
                addBamReaderInPosition(bamReader);
            }
        }

        String sequenceIntervalStr = sequenceInfo.toString();

        if(mActiveBamReaders.isEmpty())
        {
            RD_LOGGER.debug("seqRange({}) no BAM files with records found", sequenceIntervalStr, mActiveBamReaders.size());
            return;
        }

        mSamFileWriter = initialiseWriter(sequenceInfo.BamFile);

        RD_LOGGER.debug("seqRange({}) merging {} BAMs", sequenceIntervalStr, mActiveBamReaders.size());

        // begin the merge process:
        // 1. take the top record from the first BAM reader
        // 2. check the next record vs the next top record(s), and if not still top then move the BAM reader to its new position
        // 3. repeat the process
        BamSequenceReader topWriter = mActiveBamReaders.get(0);

        long recordCount = 0;

        while(!mActiveBamReaders.isEmpty())
        {
            try
            {
                mSamFileWriter.addAlignment(topWriter.current());
            }
            catch(Exception e)
            {
                RD_LOGGER.error("failed to add read from top writer({}): {}", topWriter.toString(), e.toString());
                e.printStackTrace();
                System.exit(1);
            }
            ++recordCount;

            if((recordCount % LOG_COUNT) == 0)
            {
                RD_LOGGER.trace("seqRangeId({}) merged {} records, readers(active={} finished={}) reorders({})",
                        sequenceInfo.Id, recordCount, mActiveBamReaders.size(), mFinishedBamReaders.size(), mReorderCount);
            }

            topWriter.moveNext();

            // test most likely scenario first ie that the current top BAM reader remains so
            if(topWriter.current() != null && mActiveBamReaders.size() > 1 && topWriter.isLowerOrEqualWith(mActiveBamReaders.get(1)))
            {
                continue;
            }

            mActiveBamReaders.remove(0);

            if(topWriter.finished())
            {
                RD_LOGGER.trace("seqRangeId({}) bam({}) finished", sequenceInfo.Id, topWriter.filename());
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

        RD_LOGGER.debug("seqRangeId({}) merged {} BAM files with {} records, reorder count({})",
                sequenceInfo.Id, mInputBams.size(), recordCount, mReorderCount);

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

        fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, true, new File(outputBam));
    }
}
