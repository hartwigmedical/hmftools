package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.Constants.BAM_READ_CACHE_BUFFER;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.PRIMARY;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.common.DuplicateGroup;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriter
{
    private final MarkDupsConfig mConfig;

    private final String mFilename;
    private final SAMFileWriter mBamWriter;
    private boolean mWriteSorted;

    private final SortedReadCache mReadCache;
    private final BamWriter mSharedUnsortedWriter;

    private final ReadDataWriter mReadDataWriter;

    private int mWriteCount;
    private int mNonConsensusReadCount;
    private int mConsensusReadCount;

    public BamWriter(
            final String filename, final MarkDupsConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            boolean writeSorted, final BamWriter sharedUnsortedWriter)
    {
        mFilename = filename;
        mConfig = config;
        mBamWriter = samFileWriter;
        mReadDataWriter = readDataWriter;
        mWriteSorted = writeSorted;

        if(mWriteSorted)
        {
            mReadCache = new SortedReadCache(BAM_READ_CACHE_BUFFER, (int)(config.BufferSize * 1.5), mBamWriter);
            mSharedUnsortedWriter = sharedUnsortedWriter;
        }
        else
        {
            mReadCache = null;
            mSharedUnsortedWriter = null;
        }

        mNonConsensusReadCount = 0;
        mConsensusReadCount = 0;
    }

    public String filename() { return mFilename; }

    public int recordWriteCount() { return mNonConsensusReadCount; }
    public int recordWriteCountConsensus() { return mConsensusReadCount; }
    public boolean isSorted() { return mWriteSorted; }

    public void initialiseRegion(final String chromosome, int startPosition)
    {
        if(mReadCache != null)
            mReadCache.initialiseStartPosition(chromosome, startPosition);
    }

    public void setCurrentReadPosition(int startPosition)
    {
        if(mReadCache != null)
            mReadCache.setUpperBoundPosition(startPosition);
    }

    public void onRegionComplete()
    {
        if(mReadCache != null)
            mReadCache.flush();
    }

    public synchronized void writeFragments(final List<Fragment> fragments, boolean excludeUmis)
    {
        fragments.forEach(x -> doWriteFragment(x, excludeUmis));
    }

    public synchronized void writeFragment(final Fragment fragment) { doWriteFragment(fragment, true); }

    public synchronized void writeRead(final SAMRecord read, final FragmentStatus fragmentStatus)
    {
        writeRead(read, fragmentStatus, null);
    }

    public synchronized void writeDuplicateGroup(final DuplicateGroup group, final List<SAMRecord> completeReads)
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

    private void writeRecord(final SAMRecord read)
    {
        if(mBamWriter != null)
        {
            if(mWriteSorted)
            {
                if(!mReadCache.addRecord(read))
                    mSharedUnsortedWriter.writeRecordSync(read);
            }
            else
            {
                mBamWriter.addAlignment(read);
            }
        }

        ++mWriteCount;
    }

    public synchronized void writeRecordSync(final SAMRecord read)
    {
        mBamWriter.addAlignment(read);
        ++mWriteCount;
    }

    public void close()
    {
        if(mReadCache != null)
            mReadCache.flush();

        if(mBamWriter != null)
        {
            MD_LOGGER.debug("{} records written to BAM({})", mWriteCount, mFilename);
            mBamWriter.close();
        }
    }

    @VisibleForTesting
    public void resetRecordWriteCounts()
    {
        mConsensusReadCount = 0;
        mNonConsensusReadCount = 0;
    }
}
