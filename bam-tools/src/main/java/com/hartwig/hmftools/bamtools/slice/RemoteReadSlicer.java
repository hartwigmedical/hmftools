package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class RemoteReadSlicer implements Runnable
{
    private final SliceConfig mConfig;

    private final ThreadLocal<SamReader> mSamReader;
    private final BamSlicer mBamSlicer;

    private final ReadCache mReadCache;

    private final ChrBaseRegion mCurrentSlice;
    private int mReadsProcessed = 0;

    public RemoteReadSlicer(final ChrBaseRegion slice, final SliceConfig config, final ReadCache readCache,
            final ThreadLocal<SamReader> samReader)
    {
        mConfig = config;
        mReadCache = readCache;
        mSamReader = samReader;

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mCurrentSlice = slice;
    }

    private static final int READ_LOG_COUNT = 1_000_000;

    @Override
    public void run()
    {
        BT_LOGGER.trace("remote slice({})", mCurrentSlice);
        mBamSlicer.slice(mSamReader.get(), mCurrentSlice, this::processSamRecord);
        BT_LOGGER.trace("remote slice({}) complete, processed {} reads", mCurrentSlice, mReadsProcessed);
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mReadsProcessed;

        if(mConfig.MaxRemoteReads > 0 && mReadsProcessed >= mConfig.MaxRemoteReads)
        {
            BT_LOGGER.debug("region({}) halting reads of remote region, processed {} reads",
                    mCurrentSlice, mReadsProcessed);
            mBamSlicer.haltProcessing();
            return;
        }

        if((mReadsProcessed % READ_LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("region({}) processed {} reads, current pos({})",
                    mCurrentSlice, mReadsProcessed, read.getAlignmentStart());
        }

        if(mConfig.OnlySupplementaries && !read.getSupplementaryAlignmentFlag())
            return;

        if(mConfig.DropExcluded)
        {
            // likely unmapped now with MarkDups, so not so important
            List<ChrBaseRegion> excludedRegions = ExcludedRegions.getPolyGRegions(mConfig.RefGenVersion);
            if(ChrBaseRegion.overlaps(excludedRegions, new ChrBaseRegion(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd())))
            {
                return;
            }
        }

        mReadCache.addReadRecord(read);
    }
}
