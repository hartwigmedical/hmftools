package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.io.File;
import java.util.concurrent.Callable;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RegionBamSlicer implements Callable
{
    private final SliceConfig mConfig;
    private final ChrBaseRegion mRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final ReadCache mReadCache;
    private final SliceWriter mSliceWriter;

    private int mReadsProcessed;
    private int mRemotePositionCount;

    public RegionBamSlicer(
            final ChrBaseRegion region, final SliceConfig config, final ReadCache readCache, final SliceWriter sliceWriter)
    {
        mConfig = config;
        mRegion = region;
        mReadCache = readCache;
        mSliceWriter = sliceWriter;

        mSamReader = !mConfig.RefGenomeFile.isEmpty() ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mReadsProcessed = 0;
        mRemotePositionCount = 0;
    }

    @Override
    public Long call()
    {
        BT_LOGGER.info("processing region({})", mRegion);

        mBamSlicer.slice(mSamReader, mRegion, this::processSamRecord);

        BT_LOGGER.info("region({}) complete, processed {} reads, remote positions({})",
                mRegion, mReadsProcessed, mRemotePositionCount);

        return (long)0;
    }

    private static final int LOG_COUNT = 100_000;

    @VisibleForTesting
    public void processSamRecord(final SAMRecord read)
    {
        if(!mRegion.containsPosition(read.getAlignmentStart())) // note ignores alignment end intentionally to get unmapped mates
            return;

        ++mReadsProcessed;

        if((mReadsProcessed % LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("region({}) processed {} reads, current pos({})",
                    mRegion, mReadsProcessed, read.getAlignmentStart());
        }

        if(mConfig.MaxPartitionReads > 0 && mReadsProcessed >= mConfig.MaxPartitionReads)
        {
            BT_LOGGER.debug("region({}) halting slice after {} reads", mRegion, mReadsProcessed);
            mBamSlicer.haltProcessing();
            return;
        }

        mSliceWriter.writeRead(read);

        // register any remote reads

        // check for remote mates and supplementaries
        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
        {
            checkAddRemotePosition(read.getReadName(), read.getMateReferenceName(), read.getMateAlignmentStart());
        }

        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

            checkAddRemotePosition(read.getReadName(), suppData.Chromosome, suppData.Position);
        }
    }

    private void checkAddRemotePosition(final String readId, final String chromosome, final int startPosition)
    {
        if(mRegion.containsPosition(chromosome, startPosition))
            return;

        // skip over positions already part of the initial slice
        if(mConfig.SpecificChrRegions.Regions.stream().anyMatch(x -> x.containsPosition(chromosome, startPosition)))
            return;

        mReadCache.addRemotePosition(new RemotePosition(readId, chromosome, startPosition));
        ++mRemotePositionCount;
    }
}
