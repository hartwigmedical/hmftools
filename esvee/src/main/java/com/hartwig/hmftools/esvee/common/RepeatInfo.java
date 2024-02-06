package com.hartwig.hmftools.esvee.common;

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
        String currentStr = "";
        int count = 0;
        List<RepeatInfo> repeats = null;

        int index = 0;
        while(index < bases.length - 4)
        {
            RepeatInfo repeat = findSingleBaseRepeat(bases, index);

            if(repeat == null)
            {
                repeat = findDualBaseRepeat(bases, index);
            }

            if(repeat == null)
            {
                repeat = findTripleBaseRepeat(bases, index);
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
    private static final int TRIPLE_LENGTH = 3;
    private static final int DUAL_DUAL_LENGTH = 4;

    public static RepeatInfo findSingleBaseRepeat(final byte[] bases, int index)
    {
        if(index + MIN_SINGLE_REPEAT >= bases.length)
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
        if(index + MIN_DUAL_REPEAT * 2 >= bases.length)
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
        if(index + MIN_OTHER_REPEAT * 3 >= bases.length)
            return null;

        int i = index + TRIPLE_LENGTH;
        while(i < bases.length - 2)
        {
            if(bases[i] != bases[index] || bases[i + 1] != bases[index + 1] || bases[i + 2] != bases[index + 2])
                break;

            i += TRIPLE_LENGTH;
        }

        int repeatLength = (i - index) / TRIPLE_LENGTH;

        if(repeatLength < MIN_OTHER_REPEAT)
            return null;

        String repeat = String.valueOf((char)bases[index]) + (char)bases[index + 1] + (char)bases[index + 2];
        return new RepeatInfo(index, repeat, repeatLength);
    }

    public static RepeatInfo findDualDualRepeat(final byte[] bases, int index)
    {
        if(index + MIN_OTHER_REPEAT * DUAL_DUAL_LENGTH >= bases.length)
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
}
