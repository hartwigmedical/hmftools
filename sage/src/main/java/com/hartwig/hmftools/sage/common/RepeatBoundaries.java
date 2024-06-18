package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.common.RepeatInfo.addIfUnique;
import static com.hartwig.hmftools.sage.common.RepeatInfo.extendRepeatLower;
import static com.hartwig.hmftools.sage.common.RepeatInfo.findMultiBaseRepeat;

import java.util.List;

import com.google.common.collect.Lists;

public class RepeatBoundaries
{
    public final int LowerIndex;
    public final int UpperIndex;
    public final RepeatInfo MaxRepeat;
    public final List<RepeatInfo> AllRepeats;

    public RepeatBoundaries(final int lowerIndex, final int upperIndex, final RepeatInfo maxRepeat, final List<RepeatInfo> allRepeats)
    {
        LowerIndex = lowerIndex;
        UpperIndex = upperIndex;
        MaxRepeat = maxRepeat;
        AllRepeats = allRepeats;
    }

    // set an initial search length long enough to find a min count of the longest repeat
    protected static final int REPEAT_SEARCH_LENGTH = MAX_REPEAT_LENGTH * MIN_REPEAT_COUNT;

    public static RepeatBoundaries findRepeatBoundaries(
            final byte[] bases, final int requiredIndexStart, final int requiredIndexEnd, final int maxLength, final int minCount)
    {
        int searchIndexStart = max(0, requiredIndexStart - REPEAT_SEARCH_LENGTH);

        int index = searchIndexStart;

        int minTotalRepeatLength = minCount; // being a single-base repeat

        // find the longest repeat which crosses the required indices, and a second repeat if reaches to or past the required end
        // if the longest one does not do so

        List<RepeatInfo> allRepeats = null;

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

                if(allRepeats == null)
                {
                    allRepeats = Lists.newArrayList(repeat);
                }
                else
                {
                    addIfUnique(allRepeats, repeat);
                }
            }

            ++index;
        }

        if(allRepeats == null)
            return null;

        RepeatInfo lowerRepeat = null;
        RepeatInfo upperRepeat = null;
        RepeatInfo maxRepeat = null;

        for(RepeatInfo repeat : allRepeats)
        {
            if(maxRepeat == null || repeat.totalLength() > maxRepeat.totalLength())
                maxRepeat = repeat;

            if(positionWithin(requiredIndexStart, repeat.Index, repeat.endIndex()))
            {
                if(lowerRepeat == null || repeat.Index < lowerRepeat.Index)
                    lowerRepeat = repeat;
            }

            if(positionWithin(requiredIndexEnd, repeat.Index, repeat.endIndex()))
            {
                if(upperRepeat == null || repeat.endIndex() > upperRepeat.endIndex())
                    upperRepeat = repeat;
            }
        }

        if(lowerRepeat == null)
            lowerRepeat = maxRepeat;

        if(upperRepeat == null)
            upperRepeat = maxRepeat;

        int lowerRepeatIndex = findPostRepeatIndex(lowerRepeat, bases, true);
        int upperRepeatIndex = findPostRepeatIndex(upperRepeat, bases, false);

        return new RepeatBoundaries(lowerRepeatIndex, upperRepeatIndex, maxRepeat, allRepeats);
    }

    private static int findPostRepeatIndex(final RepeatInfo repeat, final byte[] bases, boolean searchDown)
    {
        int partialBaseMatch = 0;

        if(searchDown)
        {
            int baseIndex = repeat.Index - 1;
            int indexInRepeat = repeat.repeatLength() - 1;

            while(baseIndex >= 0 && partialBaseMatch < repeat.repeatLength())
            {
                if(bases[baseIndex] != repeat.Bases.charAt(indexInRepeat))
                    break;

                ++partialBaseMatch;
                --baseIndex;
                --indexInRepeat;
            }

            return repeat.Index - partialBaseMatch - 1;
        }
        else
        {
            int baseIndex = repeat.endIndex() + 1;
            int indexInRepeat = 0;

            while(baseIndex < bases.length && partialBaseMatch < repeat.repeatLength())
            {
                if(bases[baseIndex] != repeat.Bases.charAt(indexInRepeat))
                    break;

                ++partialBaseMatch;
                ++baseIndex;
                ++indexInRepeat;
            }

            return repeat.endIndex() + partialBaseMatch + 1;
        }
    }

    public String toString()
    {
        return format("bounds(%d-%d) repeat(%s)", LowerIndex, UpperIndex, MaxRepeat != null ? MaxRepeat : "");
    }
}
