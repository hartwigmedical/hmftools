package com.hartwig.hmftools.sage.evidence;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.sage.common.ReadContext;

public class PhasedVariantGroup
{
    public final int Id;
    public final List<ReadContextCounter> PositiveReadCounters; // supported by the reads
    public final List<ReadContextCounter> NegativeReadCounters; // not supported by the reads

    public int ReadCount; // from uniquely supporting reads
    public double AllocatedReadCount; // allocated from subset groups

    // cache min and max variant positions from the positive read-counters, used for faster group matching
    private int mMinVariantPos;
    private int mMaxVariantPos;

    public PhasedVariantGroup(
            final int id, final int minVariantPos, final int maxVariantPos,
            final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters)
    {
        Id = id;
        PositiveReadCounters = posCounters;
        NegativeReadCounters = negCounters;

        mMinVariantPos = minVariantPos;
        mMaxVariantPos = maxVariantPos;

        ReadCount = 1;
        AllocatedReadCount = 0;
    }

    public int minVariantPos() { return mMinVariantPos; }
    public int maxVariantPos() { return mMaxVariantPos; }

    public boolean positionsOverlap(final PhasedVariantGroup other)
    {
        return BaseRegion.positionsOverlap(mMinVariantPos, mMaxVariantPos, other.minVariantPos(), other.maxVariantPos());
    }

    public boolean matches(final int minVariantPos, final int maxVariantPos, final List<ReadContextCounter> posCounters)
    {
        // positives need to match exactly, negatives don't
        if(minVariantPos != mMinVariantPos || maxVariantPos != mMaxVariantPos)
            return false;

        if(PositiveReadCounters.size() != posCounters.size())
            return false;

        if(PositiveReadCounters.stream().anyMatch(x -> !posCounters.contains(x)))
            return false;

        return true;
    }

    public void merge(final PhasedVariantGroup other)
    {
        ReadCount += other.ReadCount;


        // TODO: merged reads need to go into allocated not original

        AllocatedReadCount += other.AllocatedReadCount;

        int index = 0;

        // keep merged positive RCs in positional order
        for(ReadContextCounter readCounter : other.PositiveReadCounters)
        {
            boolean matched = false;

            while(index < PositiveReadCounters.size())
            {
                ReadContextCounter counter = PositiveReadCounters.get(index);
                if(counter == readCounter)
                {
                    matched = true;
                    break;
                }

                if(counter.position() <= readCounter.position())
                {
                    ++index;
                    continue;
                }
                else
                {
                    // new read counter needs to be inserted at this earlier position
                    break;
                }
            }

            if(!matched)
                PositiveReadCounters.add(index, readCounter);
        }

        // other.PositiveReadCounters.stream().filter(x -> !PositiveReadCounters.contains(x)).forEach(x -> PositiveReadCounters.add(x));
        mergeNegatives(other.NegativeReadCounters);

        mMinVariantPos = PositiveReadCounters.get(0).position();
        mMaxVariantPos = PositiveReadCounters.get(PositiveReadCounters.size() - 1).position();
    }

    public void mergeNegatives(final List<ReadContextCounter> negCounters)
    {
        negCounters.stream().filter(x -> !NegativeReadCounters.contains(x)).forEach(x -> NegativeReadCounters.add(x));
    }

    public boolean isSubsetOf(final PhasedVariantGroup other)
    {
        // returns true if this group is a subset of 'other' is a su
        if(other.PositiveReadCounters.size() < PositiveReadCounters.size())
            return false;

        if(!PositiveReadCounters.stream().allMatch(x -> other.PositiveReadCounters.contains(x)))
            return false;

        // cannot have contradictory negatives
        if(hasAnyOverlap(other.PositiveReadCounters, NegativeReadCounters) || hasAnyOverlap(PositiveReadCounters, other.NegativeReadCounters))
            return false;

        return true;
    }

    public boolean haveCommonSubset(final PhasedVariantGroup other)
    {
        // returns true if the groups contain common subsets but not all
        int countInOther = (int)PositiveReadCounters.stream().filter(x -> other.PositiveReadCounters.contains(x)).count();

        if(countInOther == 0 || countInOther == PositiveReadCounters.size() || countInOther == other.PositiveReadCounters.size())
            return false;

        // cannot have contradictory negatives
        if(hasAnyOverlap(other.PositiveReadCounters, NegativeReadCounters) || hasAnyOverlap(PositiveReadCounters, other.NegativeReadCounters))
            return false;

        return true;
    }

    private static boolean hasAnyOverlap(final List<ReadContextCounter> counters1, final List<ReadContextCounter> counters2)
    {
        return counters1.stream().anyMatch(x -> counters2.contains(x)) || counters2.stream().anyMatch(x -> counters1.contains(x));
    }

    public String toString()
    {
        return String.format("%d: range(%d - %d) pos(%d) neg(%d) rc(%d) alloc(%.1f)",
            Id, mMinVariantPos, mMaxVariantPos, PositiveReadCounters.size(), NegativeReadCounters.size(), ReadCount, AllocatedReadCount);
    }
}
