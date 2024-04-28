package com.hartwig.hmftools.common.collect;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.apache.commons.lang3.tuple.Pair;

// see https://en.wikipedia.org/wiki/Interval_tree for an explanation of this data structure
public class ImmutableIntervalTree<T>
{
    final private int mCentre;
    final private ImmutableIntervalTree<T> mLeft;
    final private ImmutableIntervalTree<T> mRight;
    final private List<Pair<BaseRegion, T>> mCentreSortedByLeft;
    final private List<Pair<BaseRegion, T>> mCentreSortedByRight;

    public ImmutableIntervalTree(final Collection<Pair<BaseRegion, T>> entries)
    {
        mCentreSortedByLeft = new ArrayList<>();
        mCentreSortedByRight = new ArrayList<>();
        if(entries.isEmpty())
        {
            mLeft = null;
            mRight = null;
            mCentre = 0;
            return;
        }

        List<Pair<BaseRegion, T>> left = new ArrayList<>();
        List<Pair<BaseRegion, T>> right = new ArrayList<>();

        Integer leftmost = null;
        Integer rightmost = null;
        for(Pair<BaseRegion, T> entry : entries)
        {
            BaseRegion interval = entry.getKey();
            leftmost = leftmost == null ? interval.start() : min(leftmost, interval.start());
            rightmost = rightmost == null ? interval.end() : max(rightmost, interval.end());
        }

        mCentre = (int) round(0.5 * ((long) leftmost + (long) rightmost));

        for(Pair<BaseRegion, T> entry : entries)
        {
            BaseRegion interval = entry.getKey();
            if(interval.end() < mCentre)
            {
                left.add(entry);
            }
            else if(interval.start() > mCentre)
            {
                right.add(entry);
            }
            else
            {
                mCentreSortedByLeft.add(entry);
                mCentreSortedByRight.add(entry);
            }
        }

        mCentreSortedByLeft.sort(Comparator.comparingInt(x -> x.getKey().start()));
        mCentreSortedByRight.sort(Comparator.comparingInt(x -> x.getKey().end()));

        mLeft = left.isEmpty() ? null : new ImmutableIntervalTree<>(left);
        mRight = right.isEmpty() ? null : new ImmutableIntervalTree<>(right);
    }

    public int size()
    {
        int output = 0;
        output += mLeft == null ? 0 : mLeft.size();
        output += mRight == null ? 0 : mRight.size();
        output += mCentreSortedByLeft.size();
        return output;
    }

    public boolean isEmpty()
    {
        if(mLeft != null && !mLeft.isEmpty())
        {
            return false;
        }

        if(mRight != null && !mRight.isEmpty())
        {
            return false;
        }

        if(!mCentreSortedByLeft.isEmpty())
        {
            return false;
        }

        return true;
    }

    private void collectContainedIntervals(final BaseRegion queryInterval, final List<Pair<BaseRegion, T>> acc)
    {
        if(queryInterval.start() < mCentre && mLeft != null && !mLeft.isEmpty())
        {
            mLeft.collectContainedIntervals(queryInterval, acc);
        }

        if(queryInterval.end() > mCentre && mRight != null && !mRight.isEmpty())
        {
            mRight.collectContainedIntervals(queryInterval, acc);
        }

        if(mCentreSortedByLeft.isEmpty())
        {
            return;
        }

        // find all centre intervals that start at the same point or after queryInterval
        Pair<BaseRegion, T> key = Pair.of(new BaseRegion(queryInterval.start(), queryInterval.start()), null);
        int leftIndex = Collections.binarySearch(mCentreSortedByLeft, key, Comparator.comparingInt(x -> x.getKey().start()));
        while(leftIndex > 0)
        {
            if(mCentreSortedByLeft.get(leftIndex - 1).getKey().start() != queryInterval.start())
            {
                break;
            }

            --leftIndex;
        }

        if(leftIndex < 0)
        {
            // index of the first element greater than the key, i.e. starting after queryInterval starts
            int insertion_point = -leftIndex - 1;
            leftIndex = insertion_point;
        }

        for(int i = leftIndex; i < mCentreSortedByLeft.size(); ++i)
        {
            Pair<BaseRegion, T> entry = mCentreSortedByLeft.get(i);
            BaseRegion centreInterval = entry.getKey();
            if(centreInterval.start() > queryInterval.end())
            {
                break;
            }

            if(centreInterval.end() <= queryInterval.end())
            {
                acc.add(entry);
            }
        }
    }

    public List<Pair<BaseRegion, T>> containedIntervals(final BaseRegion queryInterval)
    {
        List<Pair<BaseRegion, T>> acc = new ArrayList<>();
        collectContainedIntervals(queryInterval, acc);
        return acc;
    }

    public List<Pair<BaseRegion, T>> containedIntervals(int queryIntervalStart, int queryIntervalEnd)
    {
        return containedIntervals(new BaseRegion(queryIntervalStart, queryIntervalEnd));
    }
}
