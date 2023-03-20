package com.hartwig.hmftools.bamtools.slice;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class RemoteReadSlicer implements Callable
{
    private final SliceConfig mConfig;
    private final String mChromosome;
    private final List<RemotePosition> mRemotePositions;

    private final List<RemotePosition> mSlicePositions;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final SliceWriter mSliceWriter;

    private int mTotalReads;

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
    }

    public int totalReads() { return mTotalReads; }

    private static final int MAX_POSITION_DIFF = 300;

    @Override
    public Long call()
    {
        BT_LOGGER.debug("processing chromosome({}) with {} remote reads", mChromosome, mRemotePositions.size());

        for(int i = 0; i < mRemotePositions.size(); )
        {
            RemotePosition position = mRemotePositions.get(i);

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
        }

        BT_LOGGER.info("chromosome({}) remote positions({}) complete, processed {} reads",
                mChromosome, mRemotePositions.size(), mTotalReads);

        return (long)0;
    }

    private void sliceRemotePositions()
    {
        ChrBaseRegion sliceRegion = new ChrBaseRegion(
                mChromosome,
                mSlicePositions.get(0).Position,
                mSlicePositions.get(mSlicePositions.size() - 1).Position);

        mBamSlicer.slice(mSamReader, Lists.newArrayList(sliceRegion), this::processSamRecord);
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mTotalReads;

        for(int i = 0; i < mSlicePositions.size(); ++i)
        {
            RemotePosition position = mSlicePositions.get(i);

            if(position.Position == read.getAlignmentStart() && read.getReadName().equals(position.ReadId))
            {
                mSliceWriter.writeRead(read);
                mSlicePositions.remove(i);

                if(mSlicePositions.isEmpty())
                    mBamSlicer.haltProcessing();

                return;
            }
        }
    }
}
