package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.min;

import static com.hartwig.hmftools.svprep.reads.ReadFilterType.INSERT_MAP_OVERLAP;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class ReadGroup
{
    private final List<ReadRecord> mReads;
    private int mMinPositionStart;
    private boolean mIsComplete;
    private boolean mIsPartial; // all reads captured for a given region

    private short mFirstInPairExpected;
    private short mSecondInPairExpected;
    private short mFirstInPairPresent;
    private short mSecondInPairPresent;
    private boolean mHasExternalReads;

    public ReadGroup()
    {
        mReads = Lists.newArrayListWithCapacity(2);
        mFirstInPairExpected = 0;
        mSecondInPairExpected = 0;
        mFirstInPairPresent = 0;
        mSecondInPairPresent = 0;
        mHasExternalReads = false;
    }

    public final String id() { return mReads.get(0).id(); }
    public List<ReadRecord> reads() { return mReads; }

    public boolean isComplete() { return mIsComplete; }
    public boolean isPartial() { return mIsPartial; }
    public boolean isIncomplete() { return mIsPartial; }

    public boolean isJunctionFragment(boolean allowLowMapQual)
    {
        return mReads.stream().anyMatch(x -> x.filters() == 0 || (allowLowMapQual && x.filters() == ReadFilterType.MIN_MAP_QUAL.flag()));
    }

    public void addRead(final ReadRecord read, final ChrBaseRegion region)
    {
        mReads.add(read);

        if(read.isFirstOfPair())
        {
            // ++mFirstInPairExpected;
            ++mFirstInPairPresent;

            if(mFirstInPairExpected == 1)
            {
                // check whether the supplementary falls in the same region
                if(read.hasSuppAlignment())
                {
                    if(supplementaryInRegion(read.supplementaryAlignment(), region))
                        ++mFirstInPairExpected;
                    else
                        mHasExternalReads = true;
                }
            }

            if(region.containsPosition(read.MateChromosome, read.MatePosStart))
                ++mSecondInPairExpected;
            else
                mHasExternalReads = true;
        }
        else
        {
            ++mSecondInPairPresent;

            if(mSecondInPairExpected == 1)
            {
                if(read.hasSuppAlignment())
                {
                    if(supplementaryInRegion(read.supplementaryAlignment(), region))
                        ++mSecondInPairExpected;
                }
                else
                {
                    mHasExternalReads = true;
                }
            }

            if(region.containsPosition(read.MateChromosome, read.MatePosStart))
                ++mFirstInPairExpected;
            else
                mHasExternalReads = true;
        }

        if(mFirstInPairExpected == mFirstInPairPresent && mSecondInPairExpected == mSecondInPairPresent)
        {
            if(mHasExternalReads)
                mIsPartial = true;
            else
                mIsComplete = true;
        }
    }

    public String groupStatus() { return mIsComplete ? "COMPLETE" : (mIsPartial ? "PARTIAL" : "INCOMPLETE"); }

    public void merge(final ReadGroup other)
    {
        other.reads().forEach(x -> addRead(x));
    }

    private static boolean supplementaryInRegion(final SupplementaryReadData suppData, final ChrBaseRegion region)
    {
        return suppData != null && region.containsPosition(suppData.Chromosome, suppData.Position);
    }

    public ReadGroup(final ReadRecord read)
    {
        mReads = Lists.newArrayListWithCapacity(2);
        addRead(read);
    }

    public void addRead(final ReadRecord read)
    {
        mReads.add(read);

        if(mReads.size() == 1)
        {
            mMinPositionStart = read.start();
            mIsComplete = false;
            return;
        }

        mMinPositionStart = min(mMinPositionStart, read.start());
        mIsComplete = isGroupComplete(mReads);
    }

    public static boolean isGroupComplete(final List<ReadRecord> reads)
    {
        boolean firstHasSupp = false;
        boolean secondHasSupp = false;
        int firstCount = 0;
        int secondCount = 0;

        for(ReadRecord read : reads)
        {
            if(read.isFirstOfPair())
            {
                ++firstCount;
                firstHasSupp |= read.isSupplementaryAlignment() || read.supplementaryAlignment() != null;
            }
            else
            {
                ++secondCount;
                secondHasSupp |= read.isSupplementaryAlignment() || read.supplementaryAlignment() != null;
            }
        }

        return ((!firstHasSupp && firstCount == 1) || firstCount == 2) && ((!secondHasSupp && secondCount == 1) || secondCount == 2);
    }

    public int size() { return mReads.size(); }

    public int minStartPosition() { return mMinPositionStart; }

    public String toString()
    {
        return String.format("%s reads(%d) state(%s)",
                id(), mReads.size(), mIsComplete ? "complete" : (mIsPartial ? "partial" : "incomplete"));
    }
}
