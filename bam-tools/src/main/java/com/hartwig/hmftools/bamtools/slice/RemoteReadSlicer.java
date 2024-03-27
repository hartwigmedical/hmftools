package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class RemoteReadSlicer implements Callable
{
    private final SliceConfig mConfig;
    private final String mChromosome;
    private final List<RemotePosition> mRemotePositions;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final SliceWriter mSliceWriter;

    private final List<RemotePosition> mSlicePositions;
    private ChrBaseRegion mCurrentSlice;
    private int mMatchedPositions;
    private int mTotalReads;
    private int mSliceCount;

    public RemoteReadSlicer(
            final String chromosome, final List<RemotePosition> remotePositions, final SliceConfig config, final SliceWriter sliceWriter)
    {
        mConfig = config;
        mChromosome = chromosome;
        mRemotePositions = remotePositions;
        mSliceWriter = sliceWriter;

        mSamReader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mSlicePositions = Lists.newArrayList();

        mTotalReads = 0;
        mMatchedPositions = 0;
        mSliceCount = 0;
        mCurrentSlice = null;
    }

    private static final int MAX_POSITION_DIFF = 300;

    @Override
    public Long call()
    {
        BT_LOGGER.debug("processing chromosome({}) with {} remote reads", mChromosome, mRemotePositions.size());

        // likely unmapped now with MarkDups, so not so important
        List<ChrBaseRegion> excludedRegions = ExcludedRegions.getPolyGRegions(mConfig.RefGenVersion);

        for(int i = 0; i < mRemotePositions.size(); )
        {
            RemotePosition position = mRemotePositions.get(i);

            if(mConfig.DropExcluded && ChrBaseRegion.containsPosition(excludedRegions, position.Chromosome, position.Position))
            {
                ++i;
                continue;
            }

            mSlicePositions.add(position);

            int j = i + 1;

            while(j < mRemotePositions.size())
            {
                RemotePosition nextPosition = mRemotePositions.get(j);

                if(nextPosition.Position - mSlicePositions.get(mSlicePositions.size() - 1).Position > MAX_POSITION_DIFF)
                    break;

                mSlicePositions.add(nextPosition);
                ++j;
            }

            i += mSlicePositions.size();
            sliceRemotePositions();
            mSlicePositions.clear();

            if(mConfig.MaxRemoteReads > 0 && mTotalReads >= mConfig.MaxRemoteReads)
                break;
        }

        BT_LOGGER.info("chromosome({}) remote positions({}) complete, processed {} reads, slices({})",
                mChromosome, mRemotePositions.size(), mTotalReads, mSliceCount);

        return (long)0;
    }

    private void sliceRemotePositions()
    {
        mCurrentSlice = new ChrBaseRegion(
                mChromosome,
                mSlicePositions.get(0).Position,
                mSlicePositions.get(mSlicePositions.size() - 1).Position);

        BT_LOGGER.trace("remote region slice({}) for {} remote positions, matched({}/{}), processed {} reads",
                mCurrentSlice, mSlicePositions.size(), mMatchedPositions, mRemotePositions.size(), mTotalReads);

        ++mSliceCount;
        mBamSlicer.slice(mSamReader, mCurrentSlice, this::processSamRecord);
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mTotalReads;

        if(mConfig.MaxRemoteReads > 0 && mTotalReads >= mConfig.MaxRemoteReads)
        {
            BT_LOGGER.debug("chromosome({}) halting reads of remote region, matched({}/{}), processed {} reads",
                    mChromosome, mMatchedPositions, mRemotePositions.size(), mTotalReads);
            mBamSlicer.haltProcessing();
            return;
        }

        int readStartPos = read.getAlignmentStart();

        int i = 0;
        while(i < mSlicePositions.size())
        {
            RemotePosition remotePosition = mSlicePositions.get(i);

            if(readStartPos < remotePosition.Position)
                break;

            if(readStartPos > remotePosition.Position)
            {
                // suggests the read wasn't found, is unexpected
                mSlicePositions.remove(i);
                continue;
            }

            if(remotePosition.Position == readStartPos && read.getReadName().equals(remotePosition.ReadId))
            {
                mSliceWriter.writeRead(read);

                if(!expectOtherRead(read))
                    mSlicePositions.remove(i); // keep if expects a mate or supp to be in this same region

                ++mMatchedPositions;

                if(mSlicePositions.isEmpty())
                    mBamSlicer.haltProcessing();

                return;
            }

            ++i;
        }
    }

    private boolean expectOtherRead(final SAMRecord read)
    {
        if(read.getMateUnmappedFlag())
            return true;

        if(!read.getMateReferenceName().equals(mChromosome))
            return false;

        int mateStartPos = read.getMateAlignmentStart();

        return positionsOverlap(
                mCurrentSlice.start(), mCurrentSlice.end(),
                mateStartPos, mateStartPos + read.getReadBases().length * 2);
    }
}
