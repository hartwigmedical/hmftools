package com.hartwig.hmftools.cobalt.ratio;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

class DiploidRatioNormalization
{
    private int mStartIndex;
    private int mEndIndex;

    private final long mMaxWindowDistance;
    private final List<ReadRatio> mRatios;
    private final List<ReadRatio> mResults;
    private final RollingMedian mRollingMedian;

    DiploidRatioNormalization(final double expectedRatio, int maxWindowDistance, int minWindowCoverage,
            final List<ReadRatio> ratios)
    {
        mStartIndex = 0;
        mEndIndex = -1;

        mResults = Lists.newArrayList();
        mRollingMedian = new RollingMedian();

        mMaxWindowDistance = maxWindowDistance;
        mRatios = ratios;

        for(int currentIndex = 0; currentIndex < ratios.size(); currentIndex++)
        {
            final ReadRatio current = ratios.get(currentIndex);

            removeExpiredRatios(currentIndex);
            addNewRatios(currentIndex);

            double medianRatio = mRollingMedian.median();
            double correctedRatio = isValid(current) && mRollingMedian.size() >= minWindowCoverage
                    ? expectedRatio * current.ratio() / medianRatio
                    : current.ratio();

            mResults.add(ImmutableReadRatio.builder().from(current).ratio(correctedRatio).build());
        }
    }

    @NotNull
    List<ReadRatio> get()
    {
        return mResults;
    }

    private boolean isValid(@NotNull final ReadRatio ratio)
    {
        return Doubles.greaterThan(ratio.ratio(), 0);
    }

    private void addNewRatios(int currentIndex)
    {
        for(int laterIndex = mEndIndex + 1; laterIndex < mRatios.size(); laterIndex++)
        {
            final ReadRatio later = mRatios.get(laterIndex);

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

    private void addToMedian(@NotNull final ReadRatio current)
    {
        mEndIndex++;
        if(isValid(current))
        {
            mRollingMedian.add(current.ratio());
        }
    }

    private void removeExpiredRatios(int currentIndex)
    {
        for(int earlierIndex = mStartIndex; earlierIndex < currentIndex; earlierIndex++)
        {
            final ReadRatio earlier = mRatios.get(earlierIndex);
            final boolean isValid = isValid(earlier);

            if(!isValid || distance(currentIndex, earlierIndex) > mMaxWindowDistance)
            {
                if(isValid)
                {
                    mRollingMedian.remove(earlier.ratio());
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
