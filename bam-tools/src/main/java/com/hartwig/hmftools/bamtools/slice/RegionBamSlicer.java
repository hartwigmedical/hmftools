package com.hartwig.hmftools.bamtools.slice;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class RegionBamSlicer implements Runnable
{
    private final SliceParams mParams;

    final ThreadLocal<SamReader> mBamReader;
    private final BamSlicer mBamSlicer;

    private final ReadCache mReadCache;
    private final ChrBaseRegion mCurrentRegion;
    private final List<ChrBaseRegion> mLowerRegions;

    private int mReadsProcessed;

    private final TargetDepthTracker mTargetDepthTracker;

    public RegionBamSlicer(
            final ChrBaseRegion region, final List<ChrBaseRegion> allRegions, final SliceParams sliceParams, final ReadCache readCache,
            final ThreadLocal<SamReader> bamReader)
    {
        mParams = sliceParams;
        mReadCache = readCache;
        mBamReader = bamReader;

        mBamSlicer = new BamSlicer(0, !mParams.DropDuplicates, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mCurrentRegion = region;

        // make note of earlier regions to test for reads overlapping them
        mLowerRegions = allRegions.stream()
                .filter(x -> x.Chromosome.equals(mCurrentRegion.Chromosome))
                .filter(x -> x.start() < mCurrentRegion.start())
                .collect(Collectors.toList());

        mReadsProcessed = 0;

        mTargetDepthTracker = new TargetDepthTracker(
                mCurrentRegion, mParams.TargetDepth > 0, mParams.TargetDepth, mParams.ReadLength, false);
    }

    public void setTargetDepth(int targetDepth)
    {
        mTargetDepthTracker.setTargetDepth(targetDepth);
    }

    @Override
    public void run()
    {
        if(mTargetDepthTracker.applyDownsampling() && mTargetDepthTracker.targetDepth() == 0)
            return;

        mBamSlicer.slice(mBamReader.get(), mCurrentRegion, this::processSamRecord);

        // BT_LOGGER.info("region({}) complete, processed {} reads", mCurrentRegion, mReadsProcessed);
    }

    private static final int READ_LOG_COUNT = 1_000_000;

    @VisibleForTesting
    public void processSamRecord(final SAMRecord read)
    {
        if(mParams.CheckLogReads && mParams.LogReadIds.contains(read.getReadName()))
        {
            BT_LOGGER.debug("specific read({})", readToString(read));
        }

        int readStart = read.getAlignmentStart();
        int readEnd = max(read.getAlignmentEnd(), readStart); // accounting for unmapped mates

        if(!positionsOverlap(readStart, readEnd, mCurrentRegion.start(), mCurrentRegion.end()))
            return;

        // also ignore if the read overlaps with an earlier region
        if(mLowerRegions.stream().anyMatch(x -> positionsOverlap(readStart, readEnd, x.start(), x.end())))
            return;

        ++mReadsProcessed;

        if((mReadsProcessed % READ_LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("region({}) processed {} reads, current pos({})",
                    mCurrentRegion, mReadsProcessed, readStart);
        }

        if(!mTargetDepthTracker.processReadDepth(read.getDuplicateReadFlag(), readStart, readEnd))
        {
            return;
        }

        mReadCache.addReadRecord(read);

        if(mParams.MaxPartitionReads > 0 && mReadsProcessed >= mParams.MaxPartitionReads)
        {
            BT_LOGGER.debug("region({}) halting slice after {} reads", mCurrentRegion, mReadsProcessed);
            mBamSlicer.haltProcessing();
        }
    }

    @VisibleForTesting
    public int readsProcessed() { return mReadsProcessed; }




}
