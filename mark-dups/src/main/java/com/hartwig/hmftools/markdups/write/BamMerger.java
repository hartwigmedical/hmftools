package com.hartwig.hmftools.markdups.write;

import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamMerger
{
    private final String mOutputBam;
    private final List<String> mInputBams;
    private final List<BamReader> mActiveBamReaders;
    private final List<BamReader> mFinishedBamReaders;
    private final String mRefGenomeFile;
    private final SAMFileWriter mSamFileWriter;

    private int mReorderCount;

    public BamMerger(final String outputBam, final List<String> inputBams, final String refGenomeFile)
    {
        mOutputBam = outputBam;
        mInputBams = inputBams;
        mRefGenomeFile = refGenomeFile;

        mActiveBamReaders = Lists.newArrayListWithCapacity(inputBams.size());
        mFinishedBamReaders = Lists.newArrayListWithCapacity(inputBams.size());

        mSamFileWriter = initialiseWriter();
        mReorderCount = 0;
    }

    private static final int LOG_COUNT = 10_000_000;

    public boolean merge()
    {
        if(mInputBams.isEmpty() || mSamFileWriter == null)
            return false;

        // open and prepare each BAM, adding in start order
        for(String bamFile : mInputBams)
        {
            BamReader bamReader = new BamReader(bamFile);

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
            MD_LOGGER.warn("no BAM files with records found", mActiveBamReaders.size());
            return false;
        }

        MD_LOGGER.debug("merging {} BAMs", mActiveBamReaders.size());

        // begin the merge process:
        // 1. take the top record from the first BAM reader
        // 2. check the next record vs the next top record(s), and if not still top then move the BAM reader to its new position
        // 3. repeat the process
        BamReader topWriter = mActiveBamReaders.get(0);

        long recordCount = 0;

        while(!mActiveBamReaders.isEmpty())
        {
            mSamFileWriter.addAlignment(topWriter.current());
            ++recordCount;

            if((recordCount % LOG_COUNT) == 0)
            {
                MD_LOGGER.debug("merged {} records, readers(active={} finished={}) reorders({})",
                        recordCount, mActiveBamReaders.size(), mFinishedBamReaders.size(), mReorderCount);
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
                MD_LOGGER.debug("bam({}) finished", topWriter.filename());
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

        MD_LOGGER.info("merged {} BAM files with {} records, reorder count({})",
                mInputBams.size(), recordCount, mReorderCount);

        mSamFileWriter.close();

        return true;
    }

    private void addBamReaderInPosition(final BamReader bamReader)
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

    private SAMFileWriter initialiseWriter()
    {
        if(mInputBams.isEmpty() || mOutputBam == null)
            return null;

        String sampleBam = mInputBams.get(0);
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mRefGenomeFile)).open(new File(sampleBam));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        // need to check this - must be unsorted to write as
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, true, new File(mOutputBam));
    }

    private class BamReader
    {
        private final SamReader mSamReader;
        private final SAMRecordIterator mSamIterator;
        private final String mFilename;

        private SAMRecord mCurrentRecord;
        private String mCurentChromosome;
        private int mCurentChromosomeRank;
        private boolean mOnUnmappedRecords;

        public BamReader(final String bamFile)
        {
            File file = new File(bamFile);
            mFilename = file.getName();
            mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mRefGenomeFile)).open(file);
            mCurrentRecord = null;
            mCurentChromosome = "";
            mCurentChromosomeRank = -1;
            mOnUnmappedRecords = false;

            mSamIterator = mSamReader.iterator();
            moveNext();
        }

        public String filename() { return mFilename; }

        public SAMRecord current() { return mCurrentRecord; }

        public int currentPosition() { return mCurrentRecord != null ? mCurrentRecord.getAlignmentStart() : -1; }
        public int currentChromosomeRank() { return mCurentChromosomeRank; }

        public boolean isLowerOrEqualWith(final BamReader other) { return !isHigherThan(other);}

        public boolean isHigherThan(final BamReader other)
        {
            if(finished())
                return false;

            if(mOnUnmappedRecords != other.onUnmappedRecords())
                return mOnUnmappedRecords;

            if(mCurentChromosomeRank != other.currentChromosomeRank())
                return mCurentChromosomeRank > other.currentChromosomeRank();
            else
                return currentPosition() > other.currentPosition();
        }

        public SAMRecord moveNext()
        {
            if(mSamIterator.hasNext())
            {
                mCurrentRecord = mSamIterator.next();

                if(mCurrentRecord.getReadUnmappedFlag() && mCurrentRecord.getMateUnmappedFlag())
                {
                    mOnUnmappedRecords = true;
                    mCurentChromosome = "";
                    mCurentChromosomeRank = -1;
                }
                else if(!mCurentChromosome.equals(mCurrentRecord.getReferenceName()))
                {
                    mCurentChromosome = mCurrentRecord.getReferenceName();
                    mCurentChromosomeRank = HumanChromosome.chromosomeRank(mCurentChromosome);
                }
            }
            else
            {
                mCurrentRecord = null;
                mCurentChromosome = "";
                mCurentChromosomeRank = -1;
            }

            return mCurrentRecord;
        }

        public boolean finished() { return mCurrentRecord == null; }
        public boolean onUnmappedRecords() { return mOnUnmappedRecords; }

        public String toString()
        {
            String state;

            if(finished())
                state = "finished";
            else if(mOnUnmappedRecords)
                state = "unmapped";
            else
                state = format("chr(%s:%d)", mCurentChromosome, currentPosition());

            return format("%s: %s", mFilename, state);
        }
    }
}
