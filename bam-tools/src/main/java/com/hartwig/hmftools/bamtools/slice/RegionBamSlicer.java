package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bamtools.common.ReadGroup;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RegionBamSlicer implements Callable
{
    private final SliceConfig mConfig;
    private final ChrBaseRegion mRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final Map<String,ReadGroup> mReadGroupMap; // keyed by readId
    private final ReadCache mReadCache;
    private final SliceWriter mSliceWriter;

    private int mReadsProcessed;

    public RegionBamSlicer(
            final ChrBaseRegion region, final SliceConfig config, final ReadCache readCache, final SliceWriter sliceWriter)
    {
        mConfig = config;
        mRegion = region;
        mReadCache = readCache;
        mSliceWriter = sliceWriter;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mReadGroupMap = Maps.newHashMap();

        mReadsProcessed = 0;
    }

    public int totalReads() { return mReadsProcessed; }

    @Override
    public Long call()
    {
        BT_LOGGER.info("processing region({})", mRegion);

        mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);

        BT_LOGGER.info("region({}) complete, processed {} reads", mRegion, mReadsProcessed);

        return (long)0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(!positionsOverlap(mRegion.start(), mRegion.end(), read.getAlignmentStart(), read.getAlignmentEnd()))
            return;

        ++mReadsProcessed;

        if(mReadsProcessed > 0 && (mReadsProcessed % 1_000_000) == 0)
            System.gc();

        ReadGroup readGroup = mReadGroupMap.get(read.getReadName());

        if(readGroup != null)
        {
            readGroup.addRead(read);

            if(readGroup.allReadsPresent())
            {
                mReadGroupMap.remove(read.getReadName());
                mSliceWriter.writeReads(readGroup.reads());
            }
            else
            {
                List<RemotePosition> otherReadPositions = getOutstandingReadPositions(readGroup);

                if(!hasLocalRegionMatch(otherReadPositions))
                {
                    mReadGroupMap.remove(read.getReadName());
                    mSliceWriter.writeReads(readGroup.reads());

                    if(!otherReadPositions.isEmpty())
                        mReadCache.addReadGroup(otherReadPositions);
                }
            }

            return;
        }


        List<RemotePosition> otherReadPositions = getOtherReadPositions(read);

        if(hasLocalRegionMatch(otherReadPositions))
        {
            mReadGroupMap.put(read.getReadName(), new ReadGroup(read));
        }
        else
        {
            mSliceWriter.writeRead(read);
            mReadCache.addReadGroup(otherReadPositions);
        }
    }

    private boolean hasLocalRegionMatch(final List<RemotePosition> readPositions)
    {
        return readPositions.stream().anyMatch(x -> mRegion.containsPosition(x.Chromosome, x.Position));
    }

    private static List<RemotePosition> getOutstandingReadPositions(final ReadGroup readGroup)
    {
        List<RemotePosition> positions = Lists.newArrayList();

        for(SAMRecord read : readGroup.reads())
        {
            List<RemotePosition> readPositions = getOtherReadPositions(read);

            for(RemotePosition otherPosition : readPositions)
            {
                if(readGroup.reads().stream().noneMatch(x -> otherPosition.positionsMatch(x))
                && positions.stream().noneMatch(x -> x.positionsMatch(otherPosition)))
                {
                    positions.add(otherPosition);
                }
            }
        }

        return positions;
    }

    private static List<RemotePosition> getOtherReadPositions(final SAMRecord read)
    {
        List<RemotePosition> positions = org.apache.commons.compress.utils.Lists.newArrayList();

        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
        {
            positions.add(new RemotePosition(read.getReadName(), read.getMateReferenceName(), read.getMateAlignmentStart()));
        }

        SupplementaryReadData suppData = SupplementaryReadData.from(read);

        if(suppData != null)
        {
            positions.add(new RemotePosition(read.getReadName(), suppData.Chromosome, suppData.Position));
        }

        return positions;
    }
}
