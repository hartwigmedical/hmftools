package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.io.File;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RegionBamSlicer extends Thread
{
    private final SliceConfig mConfig;
    private final Queue<ChrBaseRegion> mRegionsQueue;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final ReadCache mReadCache;
    private final SliceWriter mSliceWriter;

    private ChrBaseRegion mCurrentRegion;
    private final int mTotalRegions;
    private int mReadsProcessed;
    private int mRemotePositionCount;

    public RegionBamSlicer(
            final Queue<ChrBaseRegion> regionsQueue, final SliceConfig config, final ReadCache readCache, final SliceWriter sliceWriter)
    {
        mConfig = config;
        mRegionsQueue = regionsQueue;
        mTotalRegions = regionsQueue.size();
        mReadCache = readCache;
        mSliceWriter = sliceWriter;

        mSamReader = !mConfig.RefGenomeFile.isEmpty() ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mCurrentRegion = null;

        mReadsProcessed = 0;
        mRemotePositionCount = 0;
    }

    private static final int LOG_COUNT = 1000;

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mRegionsQueue.size();
                int processedCount = mTotalRegions - remainingCount;

                mCurrentRegion = mRegionsQueue.remove();

                mBamSlicer.slice(mSamReader, mCurrentRegion, this::processSamRecord);

                BT_LOGGER.info("region({}) complete, processed {} reads, remote positions({})",
                        mCurrentRegion, mReadsProcessed, mRemotePositionCount);

                if(processedCount > 0 && (processedCount % LOG_COUNT) == 0)
                {
                    BT_LOGGER.debug("processed {} slice regions, remaining({})", processedCount, remainingCount);
                }
            }
            catch(NoSuchElementException e)
            {
                BT_LOGGER.trace("all phase tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private static final int READ_LOG_COUNT = 100_000;

    @VisibleForTesting
    public void processSamRecord(final SAMRecord read)
    {
        if(!mCurrentRegion.containsPosition(read.getAlignmentStart())) // note ignores alignment end intentionally to get unmapped mates
            return;

        ++mReadsProcessed;

        if((mReadsProcessed % READ_LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("region({}) processed {} reads, current pos({})",
                    mCurrentRegion, mReadsProcessed, read.getAlignmentStart());
        }

        if(mConfig.MaxPartitionReads > 0 && mReadsProcessed >= mConfig.MaxPartitionReads)
        {
            BT_LOGGER.debug("region({}) halting slice after {} reads", mCurrentRegion, mReadsProcessed);
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
        if(mCurrentRegion.containsPosition(chromosome, startPosition))
            return;

        // skip over positions already part of the initial slice
        if(mConfig.SpecificChrRegions.Regions.stream().anyMatch(x -> x.containsPosition(chromosome, startPosition)))
            return;

        mReadCache.addRemotePosition(new RemotePosition(readId, chromosome, startPosition));
        ++mRemotePositionCount;
    }

    @VisibleForTesting
    public void setCurrentRegion(final ChrBaseRegion region) { mCurrentRegion = region; }
}
