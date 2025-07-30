package com.hartwig.hmftools.redux.write;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.PRIMARY;

import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.redux.bqr.BaseQualRecalibration;
import com.hartwig.hmftools.redux.bqr.BaseQualityResults;
import com.hartwig.hmftools.redux.bqr.BqrRegionReader;
import com.hartwig.hmftools.redux.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.DuplicateGroupCollapser;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.common.FragmentStatus;
import com.hartwig.hmftools.redux.common.ReadInfo;

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
    private final BqrRegionReader mBqrProcessor;

    private final boolean mRecomputeFragCoords;

    protected final AtomicLong mNonConsensusReadCount;
    protected final AtomicLong mConsensusReadCount;

    public BamWriter(
            final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            @Nullable final JitterAnalyser jitterAnalyser, final BaseQualRecalibration bqr)
    {
        mFilename = filename;
        mConfig = config;
        mSamFileWriter = samFileWriter;
        mReadDataWriter = readDataWriter;
        mJitterAnalyser = jitterAnalyser;
        mBqrProcessor = new BqrRegionReader(config.Sequencing, config.RefGenome, bqr.results(), bqr.regions());

        mRecomputeFragCoords = mReadDataWriter.enabled() && (DuplicateGroupCollapser.isEnabled(mConfig.DuplicateGroupCollapse)
                || (config.Sequencing == ILLUMINA && config.UMIs.Enabled));

        mNonConsensusReadCount = new AtomicLong(0);
        mConsensusReadCount = new AtomicLong(0);
    }

    public String filename() { return mFilename; }

    public long nonConsensusWriteCount() { return mNonConsensusReadCount.get(); }
    public long consensusWriteCount() { return mConsensusReadCount.get(); }

    public long totalWriteCount() { return nonConsensusWriteCount() + consensusWriteCount(); }

    public abstract boolean isSorted();

    public void initialiseRegion(final ChrBaseRegion region)
    {
        // initialise BQR tracking
        mBqrProcessor.initialise(region);

        onRegionInitialised(region.Chromosome, region.start());
    }

    public void regionComplete()
    {
        // flush BQR results
        mBqrProcessor.onRegionComplete();

        onRegionComplete();
    }

    // methods to guide the sorted BAM writer
    public abstract void onRegionInitialised(final String chromosome, int startPosition);
    public abstract void setBoundaryPosition(int position, boolean isLower);
    public abstract void onRegionComplete();

    // the public write methods are all thread-safe, using atomic counters and then the key SAM-write methods are handled
    // in the derived sync and non-sync implementations
    public void writeNonDuplicateReads(final List<ReadInfo> readInfos)
    {
        for(ReadInfo readInfo : readInfos)
        {
            // UMIs are not captured nor written for non-duplicates
            writeRead(readInfo.read(), FragmentStatus.NONE, readInfo.coordinates().Key, "");
        }
    }

    public void writeSecondaryRead(final SAMRecord read)
    {
        writeRead(read, FragmentStatus.UNSET, "", "");
    }

    public void writeDuplicateGroup(final DuplicateGroup group)
    {
        String fragCoords = group.fragmentCoordinates().Key;
        if(group.consensusRead() != null)
        {
            SAMRecord read = group.consensusRead();
            processRecord(read);
            mConsensusReadCount.incrementAndGet();

            if(mRecomputeFragCoords)
                fragCoords = FragmentCoords.fromRead(read, false).Key;

            if(mReadDataWriter != null && mReadDataWriter.enabled())
                mReadDataWriter.writeReadData(read, PRIMARY, fragCoords, group.umiId());
        }

        // is poly-g umi collapsing the only reason this is a duplicate group?
        List<SAMRecord> remainingReads;
        if(group.readCount() - group.polyGUmiReads().size() == 1)
        {
            SAMRecord read = group.reads().get(0);
            if(mConfig.UMIs.Enabled)
                read.setAttribute(UMI_ATTRIBUTE, group.umiId());

            if(mRecomputeFragCoords)
                fragCoords = FragmentCoords.fromRead(read, false).Key;

            writeRead(read, PRIMARY, fragCoords, group.umiId());

            remainingReads = group.polyGUmiReads();
        }
        else
        {
            remainingReads = group.allReads();
        }

        for(SAMRecord read : remainingReads)
        {
            if(mConfig.UMIs.Enabled)
                read.setAttribute(UMI_ATTRIBUTE, group.umiId());

            FragmentStatus fragmentStatus = group.isPrimaryRead(read) ? PRIMARY : DUPLICATE;
            if(mRecomputeFragCoords)
                fragCoords = FragmentCoords.fromRead(read, false).Key;

            writeRead(read, fragmentStatus, fragCoords, group.umiId());
        }
    }

    protected abstract void writeRecord(final SAMRecord read);
    public abstract long unsortedWriteCount();

    protected final void processRecord(final SAMRecord read)
    {
        captureReadInfo(read);
        writeRecord(read);
    }

    public void captureReadInfo(final SAMRecord read)
    {
        if(mJitterAnalyser != null && mJitterAnalyser.bamSlicerFilter().passesFilters(read))
            mJitterAnalyser.processRead(read);

        if(mBqrProcessor.isActive())
            mBqrProcessor.processRecord(read);
    }

    public abstract void close();

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
