package com.hartwig.hmftools.markdups.write;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.BAM_READ_CACHE_BUFFER;

import com.hartwig.hmftools.markdups.MarkDupsConfig;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriterNoSync extends BamWriter
{
    private final SortedBamWriter mSortedBamWriter;
    private final BamWriter mSharedUnsortedWriter;
    private boolean mWriteSorted;

    public BamWriterNoSync(
            final String filename, final MarkDupsConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            boolean writeSorted, final BamWriter sharedUnsortedWriter)
    {
        super(filename, config, readDataWriter, samFileWriter);

        mWriteSorted = writeSorted;

        if(mWriteSorted)
        {
            int positionBuffer = (int)(config.BufferSize * 1.5);
            mSortedBamWriter = new SortedBamWriter(BAM_READ_CACHE_BUFFER, positionBuffer, samFileWriter);
            mSharedUnsortedWriter = sharedUnsortedWriter;
        }
        else
        {
            mSortedBamWriter = null;
            mSharedUnsortedWriter = null;
        }
    }

    public boolean isSorted() { return mWriteSorted; }

    public void initialiseRegion(final String chromosome, int startPosition)
    {
        if(mSortedBamWriter != null)
            mSortedBamWriter.initialiseStartPosition(chromosome, startPosition);
    }

    public void setCurrentReadPosition(int startPosition)
    {
        if(mSortedBamWriter != null)
            mSortedBamWriter.setUpperBoundPosition(startPosition);
    }

    public void onRegionComplete()
    {
        if(mSortedBamWriter != null)
            mSortedBamWriter.flush();
    }

    /*
    public void writeFragments(final List<Fragment> fragments, boolean excludeUmis)
    {
        fragments.forEach(x -> doWriteFragment(x, excludeUmis));
    }

    public void writeFragment(final Fragment fragment) { doWriteFragment(fragment, true); }

    public void writeRead(final SAMRecord read, final FragmentStatus fragmentStatus)
    {
        writeRead(read, fragmentStatus, null);
    }

    public void writeDuplicateGroup(final DuplicateGroup group, final List<SAMRecord> completeReads)
    {
        for(SAMRecord read : completeReads)
        {
            if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
            {
                writeRecord(read);
                ++mConsensusReadCount;

                mReadDataWriter.writeReadData(read, PRIMARY, group.coordinatesKey(), 0, group.id());

                continue;
            }

            if(mConfig.UMIs.Enabled)
                read.setAttribute(UMI_ATTRIBUTE, group.id());

            writeRead(read, DUPLICATE, group.coordinatesKey(), 0, group.id());
        }
    }

    private void doWriteFragment(final Fragment fragment, boolean excludeUmis)
    {
        if(excludeUmis && fragment.umi() != null) // reads in UMI groups are only written as a complete group
            return;

        if(fragment.readsWritten())
        {
            MD_LOGGER.error("fragment({}) reads already written", fragment);
            return;
        }

        fragment.setReadWritten();
        fragment.reads().forEach(x -> writeRead(x, fragment.status(), fragment));
    }

    private void writeRead(final SAMRecord read, final FragmentStatus fragmentStatus, @Nullable final Fragment fragment)
    {
        writeRead(
                read, fragmentStatus,
                fragment != null ? fragment.coordinates().Key : "",
                fragment != null ? fragment.averageBaseQual() : 0,
                fragment != null ? fragment.umi() : "");
    }

    private void writeRead(
            final SAMRecord read, final FragmentStatus fragmentStatus, final String fragmentCoordinates,
            final double avgBaseQual, final String umiId)
    {
        ++mNonConsensusReadCount;

        mReadDataWriter.writeReadData(read, fragmentStatus, fragmentCoordinates, avgBaseQual, umiId);

        read.setDuplicateReadFlag(fragmentStatus == DUPLICATE); // overwrite any existing status

        writeRecord(read);
    }
    */

    @Override
    protected void writeRecord(final SAMRecord read)
    {
        if(mSamFileWriter != null)
        {
            if(mWriteSorted)
            {
                if(mSortedBamWriter.canWriteRecord(read))
                {
                    mSortedBamWriter.addRecord(read);
                }
                else
                {
                    mSharedUnsortedWriter.writeRecord(read);
                }
            }
            else
            {
                mSamFileWriter.addAlignment(read);
            }
        }
    }

    @Override
    public void close()
    {
        if(mSortedBamWriter != null)
        {
            mSortedBamWriter.flush();

            MD_LOGGER.debug("sorted-writer records written({} writes={} avg={}) maxCache({}) to BAM({})",
                    mSortedBamWriter.written(), mSortedBamWriter.writeCount(), mSortedBamWriter.averageWriteCount(),
                    mSortedBamWriter.maxCache(), mFilename);
        }

        if(mSamFileWriter != null)
            mSamFileWriter.close();
    }
}
