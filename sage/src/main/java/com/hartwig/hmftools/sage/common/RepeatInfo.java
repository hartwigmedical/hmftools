package com.hartwig.hmftools.sage.common;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.SageConstants.MAX_REPEAT_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.common.RepeatBoundaries.REPEAT_SEARCH_LENGTH;

import java.util.List;

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

    public int totalLength() { return Count * Bases.length(); }
    public int repeatLength() { return Bases.length(); }
    public int endIndex() { return Index + totalLength() - 1; }

    public boolean matches(final RepeatInfo other)
    {
        return matches(other.Bases, other.Index, other.Count);
    }

    public boolean matches(final String bases, final int index, final int count)
    {
        return Index == index && Bases.equals(bases) && Count == count;
    }

    public boolean isMultipleOf(final RepeatInfo other)
    {
        // test if the other repeat is a longer version of this repeat
        if(other.repeatLength() < repeatLength())
            return false;

        if(other.Index < Index || other.endIndex() > endIndex())
            return false;

        for(int i = 0, j = 0; i < other.repeatLength(); ++i, ++j)
        {
            if(j >= Bases.length())
                j = 0;

            if(other.Bases.charAt(i) != Bases.charAt(j))
                return false;
        }

        return true;
    }

    public boolean isAlternativeStart(final RepeatInfo other)
    {
        // test if the other repeat is an alternative version starting later with bases in a different order
        int repeatLength = repeatLength();

        if(repeatLength == 1 || other.repeatLength() != repeatLength)
            return false;

        if(other.totalLength() > totalLength())
            return false;

        // must start between the second and last base of the other
        if(!positionWithin(other.Index, Index, Index + repeatLength - 1))
            return false;

        // check bases are an alternative cycle of the other
        for(int s = 1; s < repeatLength; ++s)
        {
            boolean allBasesMatch = true;

            for(int i = 0, j = s; i < repeatLength; ++i, ++j)
            {
                if(j >= repeatLength)
                    j = 0;

                if(other.Bases.charAt(i) != Bases.charAt(j))
                {
                    allBasesMatch = false;
                    break;
                }
            }

            if(allBasesMatch)
                return true;
        }

        return false;
    }

    public String toString() { return format("%s-%d index(%d-%d)", Bases, Count, Index, endIndex()); }

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

    public static RepeatInfo extendRepeatLower(final RepeatInfo repeatInfo, final byte[] bases)
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

            if(repeatStart - repeatLength < 0)
                break;

            repeatStart -= repeatLength;
        }

        if(extraCount == 0)
            return repeatInfo;

        return new RepeatInfo(repeatStart, repeatInfo.Bases, repeatInfo.Count + extraCount);
    }

    public static boolean addIfUnique(final List<RepeatInfo> repeats, final RepeatInfo newRepeat)
    {
        for(RepeatInfo repeat: repeats)
        {
            if(newRepeat.matches(repeat))
                return false;

            if(newRepeat.isMultipleOf(repeat) || repeat.isMultipleOf(newRepeat))
                return false;

            if(newRepeat.isAlternativeStart(repeat) || repeat.isAlternativeStart(newRepeat))
                return false;
        }

        repeats.add(newRepeat);
        return true;
    }

    public static void setReferenceMaxRepeatInfo(final SageVariant variant, final RefSequence refSequence)
    {
        int refIndex = refSequence.index(variant.position());

        int searchIndexStart = refIndex - REPEAT_SEARCH_LENGTH;
        int searchIndexEnd = refIndex + variant.ref().length();

        RepeatInfo maxRepeat = findMaxRepeat(
                refSequence.Bases, searchIndexStart, searchIndexEnd, MAX_REPEAT_LENGTH, MIN_REPEAT_COUNT, true, refIndex);

        variant.readContext().setRefMaxRepeat(maxRepeat);
    }

    public static RepeatInfo findMaxRepeat(
            final byte[] bases, final int searchIndexStart, final int searchIndexEnd,
            final int maxLength, final int minCount, boolean extendLower, final int requiredIndex) // note: required index not used in Sage
    {
        // returns the longest repeat for the defined input params based on max count not max length of the repeat sequence
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

                if(maxRepeat == null || repeat.Count > maxRepeat.Count)
                    maxRepeat = repeat;
            }

            ++index;
        }

        return maxRepeat;
    }
}
