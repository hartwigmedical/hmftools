package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.CigarUtils.getReadBoundaryPosition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
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
    private List<ChrBaseRegion> mLowerRegions;
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
        mLowerRegions = Collections.emptyList();

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

                // make note of earlier regions to test for reads overlapping them
                markLowerRegions();

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

    private void markLowerRegions()
    {
        mLowerRegions = mConfig.SliceRegions.Regions.stream()
                .filter(x -> x.Chromosome.equals(mCurrentRegion.Chromosome))
                .filter(x -> x.start() < mCurrentRegion.start())
                .collect(Collectors.toList());
    }

    private static final int READ_LOG_COUNT = 100_000;

    @VisibleForTesting
    public void processSamRecord(final SAMRecord read)
    {
        if(!positionsOverlap(read.getAlignmentStart(), read.getAlignmentEnd(), mCurrentRegion.start(), mCurrentRegion.end()))
            return;

        // also ignore if the read overlaps with an earlier region
        if(mLowerRegions.stream().anyMatch(x -> positionsOverlap(read.getAlignmentStart(), read.getAlignmentEnd(), x.start(), x.end())))
            return;

        if(mConfig.OnlySupplementaries && !read.getSupplementaryAlignmentFlag())
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
            int mateReadEnd = getMateAlignmentEnd(read);
            checkAddRemotePosition(read.getReadName(), read.getMateReferenceName(), read.getMateAlignmentStart(), mateReadEnd);
        }

        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);
            int suppPosEnd = getReadBoundaryPosition(suppData.Position, suppData.Cigar, false, false);

            checkAddRemotePosition(read.getReadName(), suppData.Chromosome, suppData.Position, suppPosEnd);
        }
    }

    private void checkAddRemotePosition(final String readId, final String chromosome, final int readPosStart, final int readPosEnd)
    {
        if(mCurrentRegion.overlaps(chromosome, readPosStart, readPosEnd))
            return;

        // skip over positions already part of the initial slice
        if(mConfig.SliceRegions.Regions.stream().anyMatch(x -> x.overlaps(chromosome, readPosStart, readPosEnd)))
            return;

        mReadCache.addRemotePosition(new RemotePosition(readId, chromosome, readPosStart));
        ++mRemotePositionCount;
    }

    @VisibleForTesting
    public void setCurrentRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;
        markLowerRegions();
    }

    @VisibleForTesting
    public int readsProcessed() { return mReadsProcessed; }
}
