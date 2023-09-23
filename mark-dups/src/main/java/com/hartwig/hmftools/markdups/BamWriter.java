package com.hartwig.hmftools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.PRIMARY;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;

import java.io.File;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.common.DuplicateGroup;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamWriter
{
    private final MarkDupsConfig mConfig;

    private final SAMFileWriter mBamWriter;
    private int mWriteCount;

    private final ReadDataWriter mReadDataWriter;

    private boolean mCacheReads;
    private final Set<String> mExpectedReads;
    private int mNonConsensusReadCount;
    private int mConsensusReadCount;

    private int mExpectedReadsSize;

    public BamWriter(final MarkDupsConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter)
    {
        mConfig = config;
        mBamWriter = samFileWriter;
        mReadDataWriter = readDataWriter;

        mCacheReads = config.runReadChecks();
        mNonConsensusReadCount = 0;
        mConsensusReadCount = 0;
        mExpectedReadsSize = 0;
        mExpectedReads = Sets.newHashSet();
    }

    public int recordWriteCount() { return mNonConsensusReadCount; }
    public int recordWriteCountConsensus() { return mConsensusReadCount; }

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

        if(mCacheReads)
            removeWrittenRead(read);
    }

    private void writeRecord(final SAMRecord read)
    {
        if(mBamWriter != null)
            mBamWriter.addAlignment(read);

        ++mWriteCount;
    }

    public void close()
    {
        MD_LOGGER.info("{} records written to BAM", mWriteCount);
        mBamWriter.close();
    }

    public synchronized void registerRead(final SAMRecord read)
    {
        if(!mCacheReads)
            return;

        mExpectedReads.add(formCachedReadString(read));
    }

    private static final int LOG_CACHE_COUNT = 100000;

    private void removeWrittenRead(final SAMRecord read)
    {
        if(!mCacheReads)
            return;

        String readStr = formCachedReadString(read);

        if(!mExpectedReads.contains(readStr))
        {
            MD_LOGGER.warn("writing uncached read({}): {}", readStr, readToString(read));
        }
        else
        {
            mExpectedReads.remove(readStr);
        }

        if(abs(mExpectedReads.size() - mExpectedReadsSize) >= LOG_CACHE_COUNT)
        {
            mExpectedReadsSize = mExpectedReads.size();
            MD_LOGGER.info("record-writer cached expected reads count({})", mExpectedReadsSize);
        }
    }

    public void logUnwrittenReads()
    {
        if(!mCacheReads || mExpectedReads.isEmpty())
            return;

        MD_LOGGER.warn("unwritten read count({})", mExpectedReads.size());

        for(String readStr : mExpectedReads)
        {
            MD_LOGGER.warn("unwritten read({})", readStr);
        }
    }

    private static String formCachedReadString(final SAMRecord read)
    {
        boolean unmapped = read.getReadUnmappedFlag();

        return format("%s_%s_%d_%s_%s_%s",
                read.getReadName(),
                unmapped ? "unmapped" : read.getReferenceName(),
                unmapped ? 0 : read.getAlignmentStart(),
                read.getReadNegativeStrandFlag() ? "fwd" : "rev",
                read.getFirstOfPairFlag() ? "R1" : "R2",
                read.getSupplementaryAlignmentFlag() ? "supp" : "prim");
    }

    @VisibleForTesting
    public void setCacheReads() { mCacheReads = true;}

    @VisibleForTesting
    public void resetRecordWriteCounts()
    {
        mConsensusReadCount = 0;
        mNonConsensusReadCount = 0;
    }
}
