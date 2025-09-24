package com.hartwig.hmftools.cobalt.ratio;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class DiploidRatioNormalization
{
    private int mStartIndex;
    private int mEndIndex;

    private final long mMaxWindowDistance;
    private final List<Double> mRatios;
    private final List<Double> mResults;
    private final RollingMedian mRollingMedian;

    public DiploidRatioNormalization(final double expectedRatio, int maxWindowDistance, int minWindowCoverage, final List<Double> ratios)
    {
        mStartIndex = 0;
        mEndIndex = -1;

        mResults = new ArrayList<>();
        mRollingMedian = new RollingMedian();

        mMaxWindowDistance = maxWindowDistance;
        mRatios = ratios;

        for(int currentIndex = 0; currentIndex < ratios.size(); currentIndex++)
        {
            final Double current = ratios.get(currentIndex);

            removeExpiredRatios(currentIndex);
            addNewRatios(currentIndex);

            double medianRatio = mRollingMedian.median();
            Double correctedRatio = current;

            if (isValid(current) && mRollingMedian.size() >= minWindowCoverage)
            {
                correctedRatio = expectedRatio * current / medianRatio;
            }

            mResults.add(correctedRatio);
        }
    }

    @NotNull
    public List<Double> get()
    {
        return mResults;
    }

    private boolean isValid(final Double ratio)
    {
        return ratio != null && Doubles.greaterThan(ratio, 0);
    }

    private void addNewRatios(final int currentIndex)
    {
        for(int laterIndex = mEndIndex + 1; laterIndex < mRatios.size(); laterIndex++)
        {
            final Double later = mRatios.get(laterIndex);

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

    private void removeExpiredRatios(final int currentIndex)
    {
        for(int earlierIndex = mStartIndex; earlierIndex < currentIndex; earlierIndex++)
        {
            final Double earlier = mRatios.get(earlierIndex);
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
