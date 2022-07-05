package com.hartwig.hmftools.svprep.reads;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.Doubles;

public class ReadRateTracker
{
    private final int mDownsampleCount;
    private final int mSegmentLength;
    private int mCurrentPosition;
    private int mReadCount;

    private boolean mIsRateLimited;
    private double mReadFactor; // how much to count each next read, defaults to 1
    private double mAdjustedReadCount;
    private int mNextReadCount;

    public ReadRateTracker(final int segmentLength, final int positionStart, final int downsampleCount)
    {
        mDownsampleCount = downsampleCount;
        mSegmentLength = segmentLength;
        mCurrentPosition = positionStart;
        mReadCount = 0;
        mReadFactor = 0;
        mAdjustedReadCount = 0;
        mIsRateLimited = false;
    }

    public int readCount() { return mReadCount; }
    public boolean isRateLimited() { return mIsRateLimited; }

    public boolean handleRead(int position)
    {
        if(position >= mCurrentPosition + mSegmentLength)
        {
            reset();
            mCurrentPosition += mSegmentLength;
        }

        ++mReadCount;

        if(!mIsRateLimited)
            return true;

        mAdjustedReadCount += mReadFactor;

        if(Doubles.greaterOrEqual(mAdjustedReadCount, mNextReadCount))
        {
            ++mNextReadCount;
            return true;
        }
        else
        {
            return false;
        }
    }

    private void reset()
    {
        if(mReadCount >= mDownsampleCount)
        {
            mIsRateLimited = true;

            double lastSegmentRate = mReadCount / (double)mDownsampleCount;
            mReadFactor = 1 / lastSegmentRate;
        }
        else
        {
            mIsRateLimited = false;
            mReadFactor = 0;
        }

        mAdjustedReadCount = 0;
        mReadCount = 0;
        mNextReadCount = 1;
    }
}
