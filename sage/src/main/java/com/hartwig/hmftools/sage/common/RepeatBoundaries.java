package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageConstants.OUTER_MIN_REPEAT_COUNT;
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

        // find all unique repeats crossing the specified core boundaries
        List<RepeatInfo> allRepeats = null;

        while(index <= min(bases.length - minTotalRepeatLength, requiredIndexEnd))
        {
            for(int repeatLength = 1; repeatLength <= maxLength; ++repeatLength)
            {
                RepeatInfo repeat = findMultiBaseRepeat(bases, index, repeatLength, minCount);

                if(repeat == null)
                    continue;

                // search backwards for longer instances of this repeat
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

        // find the longest repeat which crosses the required indices, and and secondary repeats if they reach to or past the required indices
        RepeatInfo maxRepeat = null;
        int lowerRepeatIndex = requiredIndexStart;
        int upperRepeatIndex = requiredIndexEnd;

        if(allRepeats != null)
        {
            RepeatInfo lowerRepeat = null;
            RepeatInfo upperRepeat = null;

            for(RepeatInfo repeat : allRepeats)
            {
                if(maxRepeat == null || repeat.Count > maxRepeat.Count)
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

            lowerRepeatIndex = findPostRepeatIndex(lowerRepeat, bases, true);
            upperRepeatIndex = findPostRepeatIndex(upperRepeat, bases, false);
        }

        // search for a repeat ending/starting at the outer-most base
        int outerLowerRequiredIndex = lowerRepeatIndex < requiredIndexStart ? lowerRepeatIndex : requiredIndexStart - 1;
        RepeatInfo outerLowerRepeat = findExtraLowerRepeat(bases, outerLowerRequiredIndex);

        if(outerLowerRepeat != null && outerLowerRepeat.Index < lowerRepeatIndex)
        {
            if(allRepeats == null)
                allRepeats = Lists.newArrayList();

            allRepeats.add(outerLowerRepeat);
            lowerRepeatIndex = findPostRepeatIndex(outerLowerRepeat, bases, true);
        }

        int outerUpperRequiredIndex = upperRepeatIndex > requiredIndexEnd ? upperRepeatIndex : requiredIndexEnd + 1;

        RepeatInfo outerUpperRepeat = findExtraUpperRepeat(bases, outerUpperRequiredIndex);
        if(outerUpperRepeat != null && outerUpperRepeat.endIndex() > upperRepeatIndex)
        {
            if(allRepeats == null)
                allRepeats = Lists.newArrayList();

            allRepeats.add(outerUpperRepeat);
            upperRepeatIndex = findPostRepeatIndex(outerUpperRepeat, bases, false);
        }

        return new RepeatBoundaries(lowerRepeatIndex, upperRepeatIndex, maxRepeat, allRepeats);
    }

    private static RepeatInfo findExtraLowerRepeat(final byte[] bases, final int requiredIndexStart)
    {
        // look for a novel repeat ending at the specified index and extend if found
        int searchIndexStart = max(0, requiredIndexStart - REPEAT_SEARCH_LENGTH);

        int index = searchIndexStart;
        int minTotalRepeatLength = OUTER_MIN_REPEAT_COUNT;

        // find the longest repeat which crosses the required indices, and and secondary repeats if they reach to or past the required end
        RepeatInfo maxRepeat = null;

        while(index <= min(bases.length - minTotalRepeatLength, requiredIndexStart - 1))
        {
            for(int repeatLength = 1; repeatLength <= MAX_REPEAT_LENGTH; ++repeatLength)
            {
                int requiredStart = requiredIndexStart - repeatLength * OUTER_MIN_REPEAT_COUNT;

                if(requiredStart < 0)
                    continue;

                RepeatInfo repeat = findMultiBaseRepeat(bases, index, repeatLength, MIN_REPEAT_COUNT);

                if(repeat == null)
                    continue;

                // search backwards for longer instances of this repeat
                repeat = extendRepeatLower(repeat, bases);

                if(repeat.Count < OUTER_MIN_REPEAT_COUNT)
                    continue;

                if(!positionWithin(requiredIndexStart, repeat.Index, repeat.endIndex()))
                    continue;

                if(maxRepeat == null || repeat.totalLength() > maxRepeat.totalLength())
                    maxRepeat = repeat;
            }

            ++index;
        }

        return maxRepeat;
    }

    private static RepeatInfo findExtraUpperRepeat(final byte[] bases, final int requiredIndex)
    {
        // look for a novel repeat starting at the specified index and extend if found
        int index = requiredIndex;

        RepeatInfo maxRepeat = null;

        for(int repeatLength = 1; repeatLength <= MAX_REPEAT_LENGTH; ++repeatLength)
        {
            RepeatInfo repeat = findMultiBaseRepeat(bases, index, repeatLength, OUTER_MIN_REPEAT_COUNT);

            if(repeat == null)
                continue;

            // search backwards for longer instances of this repeat
            repeat = extendRepeatLower(repeat, bases);

            if(maxRepeat == null || repeat.totalLength() > maxRepeat.totalLength())
                maxRepeat = repeat;
        }


        return maxRepeat;
    }

    private static int findPostRepeatIndex(final RepeatInfo repeat, final byte[] bases, boolean searchDown)
    {
        // extend the current repeat with any part of it found and then add 1 additional non-repeat basr
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
