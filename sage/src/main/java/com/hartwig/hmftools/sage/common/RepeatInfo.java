package com.hartwig.hmftools.sage.common;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

public class RepeatInfo
{
    public final int Index;
    public final String Bases;
    public final int Count;

    public RepeatInfo(final int index, final String bases, final int count)
    {
        Index = index;
        Bases = bases;
        Count = count;
    }

    public int length() { return Count * Bases.length(); }
    public int repeatLength() { return Bases.length(); }
    public int endIndex() { return Index + length() - 1; }

    public String toString() { return format("%d: %s-%d", Index, Bases, Count); }

    public static RepeatInfo findMaxRepeat(
            final byte[] bases, final int searchIndexStart, final int searchIndexEnd,
            final int maxLength, final int minCount, boolean extendLower, final int requiredIndex)
    {
        RepeatInfo maxRepeat = null;

        int index = searchIndexStart;

        int minTotalRepeatLength = minCount; // being a single-base repeat

        while(index <= min(bases.length - minTotalRepeatLength, searchIndexEnd))
        {
            for(int repeatLength = 1; repeatLength <= maxLength; ++repeatLength)
            {
                RepeatInfo repeat = findMultiBaseRepeat(bases, index, repeatLength, minCount);

                if(repeat == null)
                    continue;

                // search backwards for longer instances of this repeat
                if(extendLower && repeat != null)
                    repeat = extendRepeatLower(repeat, bases);

                if(requiredIndex >= 0 && !positionWithin(requiredIndex, repeat.Index, repeat.endIndex()))
                    continue;

                if(maxRepeat == null || repeat.length() > maxRepeat.length())
                    maxRepeat = repeat;
            }

            ++index;
        }

        return maxRepeat;
    }

    protected static RepeatInfo extendRepeatLower(final RepeatInfo repeatInfo, final byte[] bases)
    {
        int extraCount = 0;

        int repeatLength = repeatInfo.repeatLength();
        int repeatStart = repeatInfo.Index - repeatLength;
        byte[] repeatBases = repeatInfo.Bases.getBytes();

        while(repeatStart > 0)
        {
            boolean matched = true;
            for(int i = 0; i < repeatLength; ++i)
            {
                if(bases[repeatStart + i] != repeatBases[i])
                {
                    matched = false;
                    break;
                }
            }

            if(!matched)
            {
                if(extraCount > 0)
                    repeatStart += repeatLength; // back out last
                break;
            }

            ++extraCount;
            repeatStart -= repeatLength;
        }

        if(extraCount == 0)
            return repeatInfo;

        return new RepeatInfo(repeatStart, repeatInfo.Bases, repeatInfo.Count + extraCount);
    }

    public static RepeatInfo findMultiBaseRepeat(final byte[] bases, int index, int repeatCount, final int minCount)
    {
        if(index + (minCount * repeatCount) - 1 >= bases.length)
            return null;

        int i = index + repeatCount;
        while(i < bases.length - (repeatCount - 1))
        {
            int matchedBases = 0;

            for(int j = 0; j < repeatCount; ++j)
            {
                if(bases[i + j] == bases[index + j])
                    ++matchedBases;
                else
                    break;
            }

            if(matchedBases != repeatCount)
                break;

            i += repeatCount;
        }

        int repeatLength = (i - index) / repeatCount;

        if(repeatLength < minCount)
            return null;

        String repeat = String.valueOf((char)bases[index]);

        for(int j = 1; j < repeatCount; ++j)
            repeat += (char)bases[index + j];

        return new RepeatInfo(index, repeat, repeatLength);
    }
}
