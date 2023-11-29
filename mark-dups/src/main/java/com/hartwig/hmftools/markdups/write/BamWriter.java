package com.hartwig.hmftools.markdups.write;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.PRIMARY;

import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.markdups.MarkDupsConfig;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.common.DuplicateGroup;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public abstract class BamWriter
{
    protected final MarkDupsConfig mConfig;
    protected final String mFilename;
    protected final SAMFileWriter mSamFileWriter;
    protected final ReadDataWriter mReadDataWriter;

    protected final AtomicLong mNonConsensusReadCount;
    protected final AtomicLong mConsensusReadCount;

    public BamWriter(
            final String filename, final MarkDupsConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter)
    {
        mFilename = filename;
        mConfig = config;
        mSamFileWriter = samFileWriter;
        mReadDataWriter = readDataWriter;

        mNonConsensusReadCount = new AtomicLong(0);
        mConsensusReadCount = new AtomicLong(0);
    }

    public String filename() { return mFilename; }

    public long nonConsensusWriteCount() { return mNonConsensusReadCount.get(); }
    public long consensusWriteCount() { return mConsensusReadCount.get(); }

    public long totalWriteCount() { return nonConsensusWriteCount() + consensusWriteCount(); }

    public abstract boolean isSorted();

    // methods to guide the sorted BAM writer
    public abstract void initialiseRegion(final String chromosome, int startPosition);
    public abstract void setBoundaryPosition(int position, boolean isLower);
    public abstract void onRegionComplete();

    // the public write methods are all thread-safe, using atomic counters and then the key SAM-write methods are handled
    // in the derived sync and non-sync implementations
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
                mConsensusReadCount.incrementAndGet();

                mReadDataWriter.writeReadData(read, PRIMARY, group.coordinatesKey(), 0, group.id());

                continue;
            }

            if(mConfig.UMIs.Enabled)
                read.setAttribute(UMI_ATTRIBUTE, group.id());

            writeRead(read, DUPLICATE, group.coordinatesKey(), 0, group.id());
        }
    }

    protected abstract void writeRecord(final SAMRecord read);

    public abstract void close();

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
        mNonConsensusReadCount.incrementAndGet();

        mReadDataWriter.writeReadData(read, fragmentStatus, fragmentCoordinates, avgBaseQual, umiId);

        read.setDuplicateReadFlag(fragmentStatus == DUPLICATE); // overwrite any existing status

        writeRecord(read);
    }
}
