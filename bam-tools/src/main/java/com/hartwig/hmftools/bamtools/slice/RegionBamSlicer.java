package com.hartwig.hmftools.bamtools.slice;

import static java.lang.Math.max;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class RegionBamSlicer implements Runnable
{
    private final SliceConfig mConfig;

    final ThreadLocal<SamReader> mBamReader;
    private final BamSlicer mBamSlicer;

    private final ReadCache mReadCache;

    private final ChrBaseRegion mCurrentRegion;
    private List<ChrBaseRegion> mLowerRegions;
    private int mReadsProcessed = 0;

    public RegionBamSlicer(final ChrBaseRegion region, final SliceConfig config, final ReadCache readCache, final ThreadLocal<SamReader> bamReader)
    {
        mConfig = config;
        mReadCache = readCache;
        mBamReader = bamReader;

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mCurrentRegion = region;
        mLowerRegions = Collections.emptyList();
    }

    @Override
    public void run()
    {
        // make note of earlier regions to test for reads overlapping them
        markLowerRegions();
        mBamSlicer.slice(mBamReader.get(), mCurrentRegion, this::processSamRecord);

        // BT_LOGGER.info("region({}) complete, processed {} reads", mCurrentRegion, mReadsProcessed);
    }

    private void markLowerRegions()
    {
        mLowerRegions = mConfig.SliceRegions.Regions.stream()
                .filter(x -> x.Chromosome.equals(mCurrentRegion.Chromosome))
                .filter(x -> x.start() < mCurrentRegion.start())
                .collect(Collectors.toList());
    }

    private static final int READ_LOG_COUNT = 1_000_000;

    @VisibleForTesting
    public void processSamRecord(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();
        int readEnd = max(read.getAlignmentEnd(), readStart); // accounting for unmapped mates

        if(!positionsOverlap(readStart, readEnd, mCurrentRegion.start(), mCurrentRegion.end()))
            return;

        // also ignore if the read overlaps with an earlier region
        if(mLowerRegions.stream().anyMatch(x -> positionsOverlap(readStart, readEnd, x.start(), x.end())))
            return;

        if(mConfig.OnlySupplementaries && !read.getSupplementaryAlignmentFlag())
            return;

        ++mReadsProcessed;

        if((mReadsProcessed % READ_LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("region({}) processed {} reads, current pos({})",
                    mCurrentRegion, mReadsProcessed, readStart);
        }

        mReadCache.addReadRecord(read);

        if(mConfig.MaxPartitionReads > 0 && mReadsProcessed >= mConfig.MaxPartitionReads)
        {
            BT_LOGGER.debug("region({}) halting slice after {} reads", mCurrentRegion, mReadsProcessed);
            mBamSlicer.haltProcessing();
        }
    }

    @VisibleForTesting
    public int readsProcessed() { return mReadsProcessed; }
}
