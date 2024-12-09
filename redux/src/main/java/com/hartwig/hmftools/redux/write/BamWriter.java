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
import com.hartwig.hmftools.redux.common.ReadInfo;
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

    public BamWriter(
            final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            @Nullable final JitterAnalyser jitterAnalyser)
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
    public SAMFileWriter samFileWriter() { return mSamFileWriter; }

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
    public void writeReads(final List<ReadInfo> readInfos, boolean excludeUmis)
    {
        readInfos.forEach(x -> doWriteRead(x, excludeUmis));
    }

    public void writeRead(final SAMRecord read, final FragmentStatus fragmentStatus)
    {
        writeRead(read, fragmentStatus, null);
    }

    public void writeDuplicateGroup(final DuplicateGroup group)
    {
        if(group.consensusRead() != null)
        {
            SAMRecord read = group.consensusRead();
            processRecord(read);
            mConsensusReadCount.incrementAndGet();

            if(mReadDataWriter != null && mReadDataWriter.enabled())
                mReadDataWriter.writeReadData(read, PRIMARY, group.coordinatesKey(), group.umiId());
        }

        for(SAMRecord read : group.reads())
        {
            if(mConfig.UMIs.Enabled)
                read.setAttribute(UMI_ATTRIBUTE, group.umiId());

            FragmentStatus fragmentStatus = group.isPrimaryRead(read) ? PRIMARY : DUPLICATE;
            writeRead(read, fragmentStatus, group.coordinatesKey(), group.umiId());
        }
    }

    protected abstract void writeRecord(final SAMRecord read);
    public abstract long unsortedWriteCount();

    protected final void processRecord(final SAMRecord read)
    {
        processJitterRead(read);
        writeRecord(read);
    }

    public void processJitterRead(final SAMRecord read)
    {
        if(mJitterAnalyser != null && mJitterAnalyser.bamSlicerFilter().passesFilters(read))
        {
            mJitterAnalyser.processRead(read);
        }
    }

    public abstract void close();

    private void doWriteRead(final ReadInfo readInfo, boolean excludeDuplicates)
    {
        if(excludeDuplicates && readInfo.umi() != null) // reads in duplicate groups are only written as a complete group
            return;

        if(readInfo.readsWritten())
        {
            RD_LOGGER.error("fragment({}) reads already written", readInfo);
            return;
        }

        readInfo.setReadWritten();
        writeRead(readInfo.read(), readInfo.status(), readInfo);
    }

    private void writeRead(final SAMRecord read, final FragmentStatus fragmentStatus, @Nullable final ReadInfo readInfo)
    {
        writeRead(
                read, fragmentStatus,
                readInfo != null && readInfo.coordinates() != null ? readInfo.coordinates().Key : "",
                readInfo != null ? readInfo.umi() : "");
    }

    private void writeRead(
            final SAMRecord read, final FragmentStatus fragmentStatus, final String fragmentCoordinates, String umiId)
    {
        mNonConsensusReadCount.incrementAndGet();

        if(mReadDataWriter.enabled())
            mReadDataWriter.writeReadData(read, fragmentStatus, fragmentCoordinates, umiId);

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
