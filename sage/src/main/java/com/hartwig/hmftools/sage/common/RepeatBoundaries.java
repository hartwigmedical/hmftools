package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.common.RepeatInfo.extendRepeatLower;
import static com.hartwig.hmftools.sage.common.RepeatInfo.findMultiBaseRepeat;

public class RepeatBoundaries
{
    public final int LowerIndex;
    public final int UpperIndex;
    public final RepeatInfo MaxRepeat;

    public RepeatBoundaries(final int lowerIndex, final int upperIndex, final RepeatInfo maxRepeat)
    {
        LowerIndex = lowerIndex;
        UpperIndex = upperIndex;
        MaxRepeat = maxRepeat;
    }

    // set an initial search length long enough to find a min count of the longest repeat
    protected static final int REPEAT_SEARCH_LENGTH = MAX_REPEAT_LENGTH * MIN_REPEAT_COUNT;

    public static RepeatBoundaries findRepeatBoundaries(
            final byte[] bases, final int requiredIndexStart, final int requiredIndexEnd, final int maxLength, final int minCount)
    {
        RepeatInfo maxRepeat = null;
        RepeatInfo secondRepeat = null;

        int searchIndexStart = max(0, requiredIndexStart - REPEAT_SEARCH_LENGTH);

        int index = searchIndexStart;

        int minTotalRepeatLength = minCount; // being a single-base repeat

        // find the longest repeat which crosses the required indices, and a second repeat if reaches to or past the required end
        // if the longest one does not do so

        while(index <= min(bases.length - minTotalRepeatLength, requiredIndexEnd))
        {
            for(int repeatLength = 1; repeatLength <= maxLength; ++repeatLength)
            {
                RepeatInfo repeat = findMultiBaseRepeat(bases, index, repeatLength, minCount);

                if(repeat == null)
                    continue;

                // search backwards for longer instances of this repeat
                if(repeat != null)
                    repeat = extendRepeatLower(repeat, bases);

                if(!positionsOverlap(requiredIndexStart, requiredIndexEnd, repeat.Index, repeat.endIndex()))
                    continue;

                if(maxRepeat == null)
                {
                    maxRepeat = repeat;
                }
                else
                {
                    if(repeat.length() > maxRepeat.length())
                    {
                        maxRepeat = repeat;
                    }
                    else if(maxRepeat.endIndex() < requiredIndexEnd && repeat.endIndex() >= requiredIndexEnd)
                    {
                        if(secondRepeat == null || repeat.length() > secondRepeat.length())
                        {
                            secondRepeat = repeat;
                        }
                    }
                }
            }

            ++index;
        }

        if(maxRepeat == null)
            return null;

        int lowerRepeatIndex = maxRepeat.Index;
        int upperRepeatIndex = secondRepeat != null ? secondRepeat.endIndex() : maxRepeat.endIndex();

        return new RepeatBoundaries(lowerRepeatIndex - 1, upperRepeatIndex + 1, maxRepeat);
    }

    public String toString()
    {
        return format("bounds(%d-%d) repeat(%s)", LowerIndex, UpperIndex, MaxRepeat != null ? MaxRepeat : "");
    }
}
