package com.hartwig.hmftools.cobalt.calculations;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.cobalt.ratio.RollingMedian;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

class DiploidRatioNormaliser
{
    private int mStartIndex = 0;
    private int mEndIndex = -1;
    private Double mMedian;
    private final long mMaxWindowDistance;
    private final int mMinWindowCoverage;
    private final DescriptiveStatistics mRatioStatistics = new DescriptiveStatistics();
    private final List<BamRatio> mRatios = new ArrayList<>();
    private final RollingMedian mRollingMedian;
    private int currentIndex = 0;
    private double mExpectedRatio;

    DiploidRatioNormaliser(int maxWindowDistance, int minWindowCoverage)
    {
        mRollingMedian = new RollingMedian();
        mMaxWindowDistance = maxWindowDistance;
        mMinWindowCoverage = minWindowCoverage;
    }

    void recordRatio(final BamRatio bamRatio)
    {
        if (isValid(bamRatio.ratio()))
        {
            mRatios.add(bamRatio);
            mRatioStatistics.addValue(bamRatio.ratio());
        }
    }

    void dataCollectionFinished()
    {
        mMedian = mRatioStatistics.getPercentile(50);
    }

    void setExpectedRatio(double expectedRatio)
    {
        mExpectedRatio = expectedRatio;
    }

    double median()
    {
        checkThatDataHasBeenFinalised();
        return mMedian;
    }

    long count()
    {
        checkThatDataHasBeenFinalised();
        return mRatioStatistics.getN();
    }

    private void checkThatDataHasBeenFinalised()
    {
        if (mMedian == null)
        {
            throw new IllegalStateException("Cannot get median until data collection is finalised");
        }
    }

    double normalise(final BamRatio bamRatio)
    {
        checkThatDataHasBeenFinalised();
        if (!isValid(bamRatio.ratio()))
        {
            return bamRatio.ratio();
        }
        Preconditions.checkArgument(bamRatio.equals(mRatios.get(currentIndex)));
        final double current = bamRatio.ratio();
        removeExpiredRatios(currentIndex);
        addNewRatios(currentIndex);
        currentIndex++;

        double medianRatio = mRollingMedian.median();
        double correctedRatio = current;

        if(isValid(current) && mRollingMedian.size() >= mMinWindowCoverage)
        {
            correctedRatio = mExpectedRatio * current / medianRatio;
        }
        return correctedRatio;
    }

    private boolean isValid(final Double ratio)
    {
        return ratio != null && Doubles.greaterThan(ratio, 0);
    }

    private void addNewRatios(int currentIndex)
    {
        for(int laterIndex = mEndIndex + 1; laterIndex < mRatios.size(); laterIndex++)
        {
            final Double later = mRatios.get(laterIndex).ratio();

            if(distance(currentIndex, laterIndex) <= mMaxWindowDistance)
            {
                addToMedian(later);
            }
            else
            {
                return;
            }
        }
    }

    private void addToMedian(final Double current)
    {
        mEndIndex++;
        if(isValid(current))
        {
            mRollingMedian.add(current);
        }
    }

    private void removeExpiredRatios(int currentIndex)
    {
        for(int earlierIndex = mStartIndex; earlierIndex < currentIndex; earlierIndex++)
        {
            final double earlier = mRatios.get(earlierIndex).ratio();
            final boolean isValid = isValid(earlier);

            if(!isValid || distance(currentIndex, earlierIndex) > mMaxWindowDistance)
            {
                if(isValid)
                {
                    mRollingMedian.remove(earlier);
                }
                mStartIndex++;
            }
            else
            {
                return;
            }
        }
    }

    private long distance(int firstIndex, int secondIndex)
    {
        return Math.abs(firstIndex - secondIndex);
    }
}
