package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BamReader
{
    private final String mRefGenomeFile;
    private final List<String> mInputBams;

    private final List<BamFileReader> mBamReaders;
    private final List<BamFileReader> mActiveBamReaders;
    private final List<BamFileReader> mFinishedBamReaders;


    public BamReader(final List<String> inputBams, final String refGenomeFile)
    {
        mInputBams = inputBams;
        mRefGenomeFile = refGenomeFile;

        mBamReaders = Lists.newArrayListWithCapacity(mInputBams.size());
        mActiveBamReaders = Lists.newArrayListWithCapacity(mInputBams.size());
        mFinishedBamReaders = Lists.newArrayListWithCapacity(mInputBams.size());

        for(String bamFile : mInputBams)
        {
            mBamReaders.add(new BamFileReader(bamFile));
        }
    }

    public void sliceRegion(final ChrBaseRegion region, final Consumer<SAMRecord> consumer)
    {
        if(mBamReaders.size() == 1)
        {
            BamFileReader reader = mBamReaders.get(0);
            reader.sliceRegion(region);

            while(!reader.finished())
            {
                consumer.accept(reader.current());
                reader.moveNext();
            }

            return;
        }

        mActiveBamReaders.clear();
        mFinishedBamReaders.clear();

        for(BamFileReader reader : mBamReaders)
        {
            reader.sliceRegion(region);

            if(!reader.finished())
                addBamReaderInPosition(reader);
        }

        if(mActiveBamReaders.isEmpty())
            return;

        BamFileReader topReader = mActiveBamReaders.get(0);

        while(!mActiveBamReaders.isEmpty())
        {
            consumer.accept(topReader.current());

            topReader.moveNext();

            // test most likely scnenario first
            if(topReader.current() != null && mActiveBamReaders.size() > 1 && topReader.isLowerOrEqualWith(mActiveBamReaders.get(1)))
            {
                continue;
            }

            mActiveBamReaders.remove(0);

            if(topReader.finished())
            {
                mFinishedBamReaders.add(topReader);
            }
            else
            {
                // move the top reader to its new position in the list and repeat the process
                addBamReaderInPosition(topReader);
            }

            if(mActiveBamReaders.isEmpty())
                return;

            topReader = mActiveBamReaders.get(0);
        }
    }

    private void addBamReaderInPosition(final BamFileReader bamReader)
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

    private class BamFileReader
    {
        private final SamReader mSamReader;
        private final String mFilename;

        private SAMRecordIterator mSamIterator;
        private SAMRecord mCurrentRecord;

        public BamFileReader(final String bamFile)
        {
            File file = new File(bamFile);
            mFilename = file.getName();
            mSamReader = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.SILENT)
                    .referenceSequence(new File(mRefGenomeFile)).open(file);

            mSamIterator = null;
            mCurrentRecord = null;
        }

        public void sliceRegion(final ChrBaseRegion region)
        {
            // check that the region exists in the BAM to avoid logging an error
            if(mSamReader.getFileHeader().getSequenceDictionary().getSequences().stream()
                    .noneMatch(x -> x.getSequenceName().equals(region.Chromosome)))
            {
                return;
            }

            final QueryInterval[] queryIntervals = BamSlicer.createIntervals(List.of(region), mSamReader.getFileHeader());

            if(queryIntervals == null)
                return;

            try
            {
                mSamIterator = mSamReader.queryOverlapping(queryIntervals);
                moveNext();
            }
            catch(Exception e)
            {
                handleReadError(region.toString(), e);
            }
        }

        public void moveNext()
        {
            try
            {
                while(mSamIterator.hasNext())
                {
                    SAMRecord record = mSamIterator.next();

                    if(isMalformed(record))
                    {
                        logMalformed(record);
                        continue;
                    }

                    mCurrentRecord = record;
                    return;
                }

                mCurrentRecord = null;
                mSamIterator.close();
            }
            catch(Exception e)
            {
                handleReadError("position " + currentPosition(), e);
            }
        }

        private static boolean isMalformed(final SAMRecord record)
        {
            if(record.getReadUnmappedFlag())
                return false;

            Cigar cigar = record.getCigar();

            if(cigar.isEmpty())
                return false;

            for(CigarElement element : cigar.getCigarElements())
            {
                if(element.getLength() == 0)
                    return true;
            }

            byte[] bases = record.getReadBases();
            return bases != null && bases.length > 0 && bases.length != cigar.getReadLength();
        }

        private void logMalformed(final SAMRecord record)
        {
            RD_LOGGER.warn("BAM({}) skipping malformed read({}) cigar({}) seqLength({})",
                    mFilename, record.getReadName(), record.getCigarString(), record.getReadLength());
        }

        private void handleReadError(final String context, final Exception e)
        {
            mCurrentRecord = null;
            RD_LOGGER.warn("BAM({}) skipping unreadable read at {}: {}", mFilename, context, e.toString());
        }

        public String filename() { return mFilename; }

        public SAMRecord current() { return mCurrentRecord; }
        public int currentPosition() { return mCurrentRecord != null ? mCurrentRecord.getAlignmentStart() : -1; }
        public boolean finished() { return mCurrentRecord == null; }

        public boolean isLowerOrEqualWith(final BamFileReader other) { return !isHigherThan(other);}

        public boolean isHigherThan(final BamFileReader other)
        {
            if(finished())
                return false;

            return currentPosition() > other.currentPosition();
        }

        public String toString()
        {
            String state = finished() ? "finished" : format("position(%d)", currentPosition());
            return format("%s: %s", mFilename, state);
        }
    }
}
