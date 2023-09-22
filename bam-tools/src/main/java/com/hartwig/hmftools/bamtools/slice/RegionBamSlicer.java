package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.io.File;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
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

        mSamReader = !mConfig.RefGenomeFile.isEmpty() ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

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

    private static final int LOG_COUNT = 100_000;

    @VisibleForTesting
    public void processSamRecord(final SAMRecord read)
    {
        if(!positionsOverlap(mRegion.start(), mRegion.end(), read.getAlignmentStart(), read.getAlignmentEnd()))
            return;

        ++mReadsProcessed;

        if((mReadsProcessed % LOG_COUNT) == 0)
        {
            BT_LOGGER.debug("region({}) processed {} reads, current pos({}) readGroupMap({})",
                    mRegion, mReadsProcessed, read.getAlignmentStart(), mReadGroupMap.size());
        }

        if(mConfig.MaxPartitionReads > 0 && mReadsProcessed >= mConfig.MaxPartitionReads)
        {
            BT_LOGGER.debug("region({}) halting slice after {} reads", mRegion, mReadsProcessed);
            mBamSlicer.haltProcessing();
            return;
        }

        mSliceWriter.writeRead(read);

        // register any remote reads

        ReadGroup readGroup = mReadGroupMap.get(read.getReadName());

        boolean groupExists = readGroup != null;

        // check for remote mates and supplementaries
        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
        {
            if(mRegion.containsPosition(read.getMateReferenceName(), read.getMateAlignmentStart()))
            {
                if(readGroup == null)
                    readGroup = new ReadGroup(read.getReadName());

                readGroup.registerMate();
            }
            else
            {
                mReadCache.addRemotePosition(new RemotePosition(read.getReadName(), read.getMateReferenceName(), read.getMateAlignmentStart()));
            }
        }

        if(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
        {
            SupplementaryReadData suppData = SupplementaryReadData.from(read);

            if(mRegion.containsPosition(suppData.Chromosome, suppData.Position))
            {
                if(readGroup == null)
                    readGroup = new ReadGroup(read.getReadName());

                readGroup.registerSupplementaryData(read, suppData);
            }
            else
            {
                mReadCache.addRemotePosition(new RemotePosition(read.getReadName(), suppData.Chromosome, suppData.Position));
            }
        }

        if(readGroup == null)
            return;

        readGroup.registerRead(read);

        if(!groupExists && !readGroup.localReadsReceived())
        {
            mReadGroupMap.put(read.getReadName(), readGroup);
        }
        else if(groupExists && readGroup.localReadsReceived())
        {
            mReadGroupMap.remove(read.getReadName());
        }
    }

    @VisibleForTesting
    public Map<String,ReadGroup> readGroupMap() { return mReadGroupMap; }

    /*
    private boolean hasLocalRegionMatch(final List<RemotePosition> readPositions)
    {
        return readPositions.stream().anyMatch(x -> mRegion.containsPosition(x.Chromosome, x.Position));
    }

    private List<RemotePosition> getRemoteReadPositions(final ReadGroup readGroup)
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

    private List<RemotePosition> getOtherReadPositions(final SAMRecord read)
    {
        List<RemotePosition> positions = Lists.newArrayList();

        if(read.getReadPairedFlag() && !read.getMateUnmappedFlag())
        {
            positions.add(new RemotePosition(read.getReadName(), read.getMateReferenceName(), read.getMateAlignmentStart()));
        }

        if(!read.getSupplementaryAlignmentFlag() && mConfig.DropRemoteSupplementaries)
            return positions;

        SupplementaryReadData suppData = SupplementaryReadData.from(read);

        if(suppData != null)
        {
            positions.add(new RemotePosition(read.getReadName(), suppData.Chromosome, suppData.Position));
        }

        return positions;
    }
    */
}
