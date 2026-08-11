package com.hartwig.hmftools.bamtools.slice;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class TargetDepthTracker
{
    // down-sampling params and state
    private final ChrBaseRegion mRegion;
    private final boolean mApplyDownsampling;
    private final int mReadLength;
    private final boolean mUseCoveragePerWindow;

    private int mWindowSize;
    private int mTargetDepth;

    private int mWindowStartPosition;
    private int mWindowReadCount; // count of reads starting in the current window
    private int mMaxReadsPerWindow;
    private int mBaseMaxReadsPerWindow;

    private final List<DepthWindow> mDepthWindows;

    private static final double DEPTH_WINDOW_FRACTION = 0.25;
    private static final double ROLLING_WINDOW_FRACTION = 0.5;
    private static final double MAX_CARRY_OVER_PERCENT = 0.25;

    public TargetDepthTracker(
            final ChrBaseRegion region, boolean applyDownsampling, int targetDepth, int readLength,  boolean useCoveragePerWindow)
    {
        mRegion = region;
        mApplyDownsampling = applyDownsampling;
        mReadLength = readLength;
        mUseCoveragePerWindow = useCoveragePerWindow;

        mWindowStartPosition = 0;
        mMaxReadsPerWindow = 0;
        mWindowReadCount = 0;

        mDepthWindows = Lists.newArrayList();
        mWindowSize = 0;

        setTargetDepth(targetDepth);
    }

    public void setTargetDepth(int targetDepth)
    {
        mTargetDepth = targetDepth;

        if(mUseCoveragePerWindow)
            computeDownsamplingFactors();
        else
            computeWindowReadCounts();
    }

    private boolean isSmallRegion() { return mReadLength >= mRegion.baseLength(); }

    public int targetDepth() { return mTargetDepth; }
    public boolean applyDownsampling() { return mApplyDownsampling; }

    public int windowSize() { return mWindowSize; }
    public int maxReadsPerWindow() { return mMaxReadsPerWindow; }
    public int baseMaxReadsPerWindow() { return mBaseMaxReadsPerWindow; }
    public int windowReadCount() { return mWindowReadCount; }

    private void computeWindowReadCounts()
    {
        mWindowReadCount = 0;
        mWindowStartPosition = mRegion.start();

        if(isSmallRegion())
        {
            // no need for a rolling window for very short regions
            mWindowSize = mRegion.baseLength();
            mMaxReadsPerWindow = mBaseMaxReadsPerWindow = mTargetDepth;
            return;
        }

        mWindowSize = (int)ceil(mReadLength * ROLLING_WINDOW_FRACTION);
        mBaseMaxReadsPerWindow = (int)ceil(mTargetDepth * (double)mWindowSize / mReadLength);
        mMaxReadsPerWindow = mBaseMaxReadsPerWindow;
    }

    public void computeDownsamplingFactors()
    {
        if(!mApplyDownsampling)
            return;

        mDepthWindows.clear();
        mWindowSize = (int)round(mReadLength * DEPTH_WINDOW_FRACTION);

        int minLength = min(mReadLength, mRegion.baseLength());
        int requiredWindows = (int)ceil(minLength / mWindowSize) + 1;

        int windowStart = mRegion.start();
        for(int i = 0; i < requiredWindows; ++i)
        {
            int windowEnd = windowStart + mWindowSize - 1;
            mDepthWindows.add(new DepthWindow(windowStart, windowEnd));
            windowStart = windowEnd + 1;
        }
    }

    public boolean processReadDepth(boolean isDuplicate, int readStart, int readEnd)
    {
        if(!mApplyDownsampling)
            return true;

        if(readStart > mRegion.end() || readEnd < mRegion.start())
            return true;

        if(mUseCoveragePerWindow)
            return processByWindowCoverage(isDuplicate, readStart, readEnd);
        else
            return processByRollingWindow(isDuplicate, readStart);
    }

    private boolean processByRollingWindow(boolean isDuplicate, int readStart)
    {
        if(!isSmallRegion())
        {
            int carryOverReads = 0;
            int windowEnd = mWindowStartPosition + mWindowSize - 1;

            if(readStart > windowEnd)
            {
                if(readStart < windowEnd + mWindowSize)
                {
                    int spareCapacity = mBaseMaxReadsPerWindow - mWindowReadCount;
                    carryOverReads = min(spareCapacity, (int)round(MAX_CARRY_OVER_PERCENT * mBaseMaxReadsPerWindow));
                }

                mWindowStartPosition = readStart;
                mWindowReadCount = 0;
                mMaxReadsPerWindow = mBaseMaxReadsPerWindow + carryOverReads;
            }
        }

        if(mWindowReadCount >= mMaxReadsPerWindow)
            return false;

        if(!isDuplicate)
            ++mWindowReadCount;

        return true;
    }

    private boolean processByWindowCoverage(boolean isDuplicate, int readStart, int readEnd)
    {
        // purge windows once read start has progressed past their end
        while(!mDepthWindows.isEmpty() && readStart > mDepthWindows.get(0).PosEnd)
        {
            mDepthWindows.remove(0);
        }

        // add new windows if required
        int maxPosEnd = min(readEnd, mRegion.end());
        int lastWindowIndex = mDepthWindows.size() - 1;
        if(mDepthWindows.isEmpty() || maxPosEnd > mDepthWindows.get(lastWindowIndex).PosEnd)
        {
            int windowEnd = !mDepthWindows.isEmpty() ? mDepthWindows.get(lastWindowIndex).PosEnd : 0;
            int windowStart = !mDepthWindows.isEmpty() ? mDepthWindows.get(lastWindowIndex).PosEnd + 1 : readStart;

            while(maxPosEnd > windowEnd)
            {
                windowEnd = windowStart + mWindowSize - 1;
                mDepthWindows.add(new DepthWindow(windowStart, windowEnd));
                windowStart = windowEnd + 1;
            }
        }

        boolean anyAlreadyFull = false;

        for(DepthWindow depthWindow : mDepthWindows)
        {
            if(positionsOverlap(readStart, readEnd, depthWindow.PosStart, depthWindow.PosEnd))
            {
                // duplicates are included but not counted towards down-sampling so the intended depth is based on primaries
                anyAlreadyFull |= depthWindow.isFull();

                if(!depthWindow.isFull() && !isDuplicate)
                    depthWindow.addRead(readStart, readEnd);
            }
        }

        return !anyAlreadyFull;
    }

    private class DepthWindow
    {
        public final int PosStart;
        public final int PosEnd;

        private int mReadCount;
        private int mTotalCoveredBases;
        private boolean mFull;

        public DepthWindow(final int posStart, final int posEnd)
        {
            PosStart = posStart;
            PosEnd = posEnd;
            mReadCount = 0;
            mTotalCoveredBases = 0;
            mFull = false;
        }

        public void addRead(int alignmentStart, int alignentEnd)
        {
            int overlap = min(PosEnd, alignentEnd) - max(PosStart, alignmentStart) + 1;

            if(overlap > 0)
            {
                ++mReadCount;
                mTotalCoveredBases += overlap;

                mFull = averageDepth() >= mTargetDepth;
            }
        }

        public boolean isFull() { return mFull; }
        public double averageDepth() { return (int)(mTotalCoveredBases / (double)(PosEnd - PosStart + 1)); }

        public String toString()
        {
            return format("window(%d - %d) reads(%d) avgDepth(%.1f)", PosStart, PosEnd, mReadCount, averageDepth());
        }
    }

}
