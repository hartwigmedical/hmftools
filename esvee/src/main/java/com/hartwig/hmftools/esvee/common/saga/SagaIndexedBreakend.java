package com.hartwig.hmftools.esvee.common.saga;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.utils.IntPair;

// For matching to SAGA variants by genomic location.
public record SagaIndexedBreakend(
        SagaBreakend breakend,
        SagaVariant variant
) implements Comparable<SagaIndexedBreakend>
{
    public String chromosome()
    {
        return breakend.chromosome();
    }

    public int position()
    {
        return breakend.position();
    }

    public Orientation orientation()
    {
        return breakend.orientation();
    }

    @Override
    public int compareTo(final SagaIndexedBreakend other)
    {
        return Integer.compare(position(), other.position());
    }

    // Finds the position in the sorted list of breakends for a single chromosome.
    // Returns a range of [startIndex, endIndex) of breakends with the specified position.
    // If there are no breakends with the position, returns [index, index) where index is the insertion position.
    static IntPair binarySearch(final List<SagaIndexedBreakend> breakends, int position)
    {
        // Only the position is compared, so the other values are simply placeholders.
        SagaIndexedBreakend target = new SagaIndexedBreakend(new SagaBreakend(new BasePosition("", position), Orientation.FORWARD), null);
        int index = Collections.binarySearch(breakends, target);
        if(index >= 0)
        {
            // Exact position found.

            // Need to search backwards to the first breakend with the same position.
            int startIndex = index;
            while(true)
            {
                int prevIndex = startIndex - 1;
                if(prevIndex >= 0 && breakends.get(prevIndex).position() == position)
                {
                    startIndex = prevIndex;
                }
                else
                {
                    break;
                }
            }

            // Need to search forwards to past the last breakend with the same position.
            int endIndex = index;
            while(true)
            {
                int nextIndex = endIndex + 1;
                if(nextIndex < breakends.size() && breakends.get(nextIndex).position() == position)
                {
                    endIndex = nextIndex;
                }
                else
                {
                    break;
                }
            }
            // Move to one past the end to form a half-open range.
            ++endIndex;

            return new IntPair(startIndex, endIndex);
        }
        else
        {
            // The exact position was not found.
            // This is the index of the breakend with the next greater position.
            index = -(index + 1);
            return new IntPair(index, index);
        }
    }
}
