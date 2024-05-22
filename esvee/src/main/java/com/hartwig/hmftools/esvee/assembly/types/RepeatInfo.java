package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

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

    public int postRepeatIndex() { return Index + length(); }
    public int length() { return Count * Bases.length(); }

    public boolean matchesType(final RepeatInfo other) { return Bases.equals(other.Bases); }

    public String toString() { return format("%d: %s-%d", Index, Bases, Count); }

    public static final int MIN_SINGLE_REPEAT = 4;
    public static final int MIN_DUAL_REPEAT = 3;
    public static final int MIN_OTHER_REPEAT = 2;

    public static List<RepeatInfo> findRepeats(final byte[] bases)
    {
        // types of repeats, single, dual, triples (ATCATC) and dual x2 (AATTAATT)
        // favour longer repeat types
        List<RepeatInfo> repeats = null;

        int index = 0;
        while(index <= bases.length - MIN_SINGLE_REPEAT)
        {
            RepeatInfo repeat = findSingleBaseRepeat(bases, index);

            if(repeat == null)
            {
                repeat = findDualBaseRepeat(bases, index);
            }

            if(repeat == null)
            {
                repeat = findMultiBaseRepeat(bases, index, THREE_LENGTH);
            }

            if(repeat == null)
            {
                repeat = findMultiBaseRepeat(bases, index, FOUR_LENGTH);
            }

            if(repeat == null)
            {
                repeat = findMultiBaseRepeat(bases, index, FIVE_LENGTH);
            }

            if(repeat == null)
            {
                repeat = findDualDualRepeat(bases, index);
            }

            if(repeat != null)
            {
                if(repeats == null)
                    repeats = Lists.newArrayList(repeat);
                else
                    repeats.add(repeat);

                index += repeat.length();
            }
            else
            {
                ++index;
            }
        }

        return repeats;
    }

    public static String repeatsAsString(final List<RepeatInfo> repeats)
    {
        if(repeats.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(" ");
        repeats.forEach(x -> sj.add(format("%dx%s", x.Count, x.Bases)));

        return format("%d %s", repeats.size(), sj);
    }

    private static final int DUAL_LENGTH = 2;
    private static final int THREE_LENGTH = 3;
    private static final int FOUR_LENGTH = 4;
    private static final int FIVE_LENGTH = 5;
    private static final int DUAL_DUAL_LENGTH = 4;

    public static RepeatInfo findSingleBaseRepeat(final byte[] bases, int index)
    {
        // the first base counts towards the repeat, so if length = 10, last index = 9, index cannot be higher than 9 - 4 + 1 = 6
        if(index + MIN_SINGLE_REPEAT - 1 >= bases.length)
            return null;

        int i = index + 1;
        while(i < bases.length)
        {
            if(bases[i] != bases[index])
                break;

            ++i;
        }

        int repeatLength = i - index;

        return repeatLength >= MIN_SINGLE_REPEAT ? new RepeatInfo(index, String.valueOf((char)bases[index]), repeatLength) : null;
    }

    public static RepeatInfo findDualBaseRepeat(final byte[] bases, int index)
    {
        if(index + MIN_DUAL_REPEAT * 2 - 1 >= bases.length)
            return null;

        int i = index + DUAL_LENGTH;
        while(i < bases.length - 1)
        {
            if(bases[i] != bases[index] || bases[i + 1] != bases[index + 1])
                break;

            i += DUAL_LENGTH;
        }

        int repeatLength = (i - index) / DUAL_LENGTH;

        if(repeatLength < MIN_DUAL_REPEAT)
            return null;

        String repeat = String.valueOf((char)bases[index]) + (char)bases[index + 1];
        return new RepeatInfo(index, repeat, repeatLength);
    }

    public static RepeatInfo findTripleBaseRepeat(final byte[] bases, int index)
    {
        return findMultiBaseRepeat(bases, index, THREE_LENGTH);
    }

    public static RepeatInfo findMultiBaseRepeat(final byte[] bases, int index, int repeatCount)
    {
        if(index + MIN_OTHER_REPEAT * repeatCount - 1 >= bases.length)
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

        if(repeatLength < MIN_OTHER_REPEAT)
            return null;

        String repeat = String.valueOf((char)bases[index]);

        for(int j = 1; j < repeatCount; ++j)
            repeat += (char)bases[index + j];

        return new RepeatInfo(index, repeat, repeatLength);
    }

    public static RepeatInfo findDualDualRepeat(final byte[] bases, int index)
    {
        if(index + MIN_OTHER_REPEAT * DUAL_DUAL_LENGTH - 1 >= bases.length)
            return null;

        // check for initial repeat type
        if(bases[index] != bases[index + 1] || bases[index + 2] != bases[index + 3])
            return null;

        int i = index + DUAL_DUAL_LENGTH;
        while(i < bases.length - 3)
        {
            if(bases[i] != bases[index] || bases[i + 1] != bases[index + 1] || bases[i + 2] != bases[index + 2]
            || bases[i + 3] != bases[index + 3])
            {
                break;
            }

            i += 4;
        }

        int repeatLength = (i - index) / DUAL_DUAL_LENGTH;

        if(repeatLength < MIN_OTHER_REPEAT)
            return null;

        String repeat = String.valueOf((char)bases[index]) + (char)bases[index + 1] + (char)bases[index + 2] + (char)bases[index + 3];
        return new RepeatInfo(index, repeat, repeatLength);
    }

    public static String buildTrimmedRefBaseSequence(final JunctionAssembly assembly, final int maxSequenceLength)
    {
        if(assembly.repeatInfo().isEmpty())
            return "";

        int refBaseLength = assembly.refBaseLength();
        int seqStart = assembly.isForwardJunction() ? assembly.junctionIndex() - refBaseLength + 1 : assembly.junctionIndex();
        int seqEnd = assembly.isForwardJunction() ? assembly.junctionIndex() : assembly.junctionIndex() + refBaseLength - 1;

        StringBuilder refBasesTrimmed = new StringBuilder();
        int currentIndex = seqStart;

        List<RepeatInfo> repeats = assembly.repeatInfo();
        int currentRepeatIndex = 0;
        int trimmedBasesLength = 0;
        RepeatInfo currentRepeat = repeats.get(currentRepeatIndex);

        while(currentIndex <= seqEnd)
        {
            while(currentRepeat != null && currentRepeat.Index < currentIndex)
            {
                ++currentRepeatIndex;

                if(currentRepeatIndex >= repeats.size())
                {
                    currentRepeat = null;
                    break;
                }
                else
                {
                    currentRepeat = repeats.get(currentRepeatIndex);
                }
            }

            if(currentRepeat != null && currentRepeat.Index == currentIndex)
            {
                if(trimmedBasesLength != 0)
                    refBasesTrimmed.append("_");

                refBasesTrimmed.append(format("%s%d_", currentRepeat.Bases, currentRepeat.Count));
                currentIndex += currentRepeat.length();
                trimmedBasesLength += currentRepeat.Bases.length() * 2;
            }
            else
            {
                refBasesTrimmed.append((char)assembly.bases()[currentIndex]);
                ++currentIndex;
                ++trimmedBasesLength;
            }

            if(trimmedBasesLength >= maxSequenceLength)
                break;
        }

        return refBasesTrimmed.toString();
    }

    public static int calcTrimmedBaseLength(final int seqStart, final int seqEnd, final List<RepeatInfo> repeats)
    {
        if(repeats == null || repeats.isEmpty())
            return seqEnd - seqStart + 1;

        int currentIndex = seqStart;

        int currentRepeatIndex = 0;
        int trimmedBasesLength = 0;
        RepeatInfo currentRepeat = repeats.get(currentRepeatIndex);

        while(currentIndex <= seqEnd)
        {
            while(currentRepeat != null && currentRepeat.Index < currentIndex)
            {
                ++currentRepeatIndex;

                if(currentRepeatIndex >= repeats.size())
                {
                    currentRepeat = null;
                    break;
                }
                else
                {
                    currentRepeat = repeats.get(currentRepeatIndex);
                }
            }

            if(currentRepeat != null && currentRepeat.Index == currentIndex)
            {
                currentIndex += currentRepeat.length();
                trimmedBasesLength += currentRepeat.Bases.length() * 2; // only count the repeat twice
            }
            else
            {
                ++currentIndex;
                ++trimmedBasesLength;
            }
        }

        return trimmedBasesLength;
    }
}
