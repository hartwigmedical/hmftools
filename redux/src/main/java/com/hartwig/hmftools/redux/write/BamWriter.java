package com.hartwig.hmftools.redux.write;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.PRIMARY;

import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.Fragment;
import com.hartwig.hmftools.redux.common.FragmentStatus;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public abstract class BamWriter
{
    protected final ReduxConfig mConfig;
    protected final String mFilename;
    protected final SAMFileWriter mSamFileWriter;
    protected final ReadDataWriter mReadDataWriter;
    private final JitterAnalyser mJitterAnalyser;

    protected final AtomicLong mNonConsensusReadCount;
    protected final AtomicLong mConsensusReadCount;

    public BamWriter(final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter,
            final SAMFileWriter samFileWriter, @Nullable final JitterAnalyser jitterAnalyser)
    {
        mFilename = filename;
        mConfig = config;
        mSamFileWriter = samFileWriter;
        mReadDataWriter = readDataWriter;
        mJitterAnalyser = jitterAnalyser;

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
                processRecord(read);
                mConsensusReadCount.incrementAndGet();

                mReadDataWriter.writeReadData(read, PRIMARY, group.coordinatesKey(), 0, group.umiId());

                continue;
            }

            if(mConfig.UMIs.Enabled)
                read.setAttribute(UMI_ATTRIBUTE, group.umiId());

            writeRead(read, DUPLICATE, group.coordinatesKey(), 0, group.umiId());
        }
    }

    protected abstract void writeRecord(final SAMRecord read);

    protected final void processRecord(final SAMRecord read)
    {
        if(mJitterAnalyser != null && mJitterAnalyser.bamSlicerFilter().passesFilters(read))
        {
            mJitterAnalyser.processRead(read);
        }

        writeRecord(read);
    }

    public abstract void close();

    private void doWriteFragment(final Fragment fragment, boolean excludeDuplicates)
    {
        if(excludeDuplicates && fragment.umi() != null) // reads in duplicate groups are only written as a complete group
            return;

        if(fragment.readsWritten())
        {
            RD_LOGGER.error("fragment({}) reads already written", fragment);
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

        if(fragmentStatus == DUPLICATE)
        {
            if(mConfig.DropDuplicates)
                return;

            read.setDuplicateReadFlag(true); // overwrite any existing status
        }

        processRecord(read);
    }

    public String toString()
    {
        return format("file(%s)", FileWriterUtils.filenamePart(mFilename));
    }
}
