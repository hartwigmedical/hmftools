package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REPEAT_MIN_COUNT;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class RepeatInfo
{
    public final int Index;
    public final String Bases;
    public final int Count;

    private final Orientation mDirection; // direction found: if -ve then Index is at the end, and lastIndex before Index

    public RepeatInfo(final int index, final String bases, final int count)
    {
        this(index, bases, count, null); // assumed forward
    }

    public RepeatInfo(final int index, final String bases, final int count, final Orientation orientation)
    {
        Index = index;
        Bases = bases;
        Count = count;
        mDirection = orientation;
    }

    public int lastIndex()
    {
        return mDirection == null || mDirection.isForward() ? indexEnd() : indexStart();
    }

    public int postRepeatIndex()
    {
        return mDirection == null || mDirection.isForward() ? Index + totalLength() : Index - totalLength();
    }

    public int indexStart() { return direction().isForward() ? Index : Index - totalLength() + 1; }
    public int indexEnd() { return direction().isForward() ? Index + totalLength() - 1 : Index; }

    public int totalLength() { return Count * Bases.length(); }
    public int repeatLength() { return Bases.length(); }
    public Orientation direction() { return mDirection != null ? mDirection : Orientation.FORWARD; }

    public boolean matchesType(final RepeatInfo other) { return Bases.equals(other.Bases); }

    public boolean overlaps(final RepeatInfo other)
    {
        return positionsOverlap(indexStart(), indexEnd(), other.indexStart(), other.indexEnd());
    }

    public String toString() { return format("%d-%d: %s-%d", indexStart(), indexEnd(), Bases, Count); }

    public static List<RepeatInfo> findRepeats(final byte[] bases)
    {
        return findRepeats(bases, REPEAT_MIN_COUNT);
    }

    public static List<RepeatInfo> findRepeats(final byte[] bases, final int minRepeatCount)
    {
        // types of repeats, single, dual, triples (ATCATC) and dual x2 (AATTAATT)
        // favour longer repeat types
        List<RepeatInfo> repeats = null;

        int index = 0;
        while(index <= bases.length - REPEAT_MIN_COUNT)
        {
            RepeatInfo repeat = findSingleBaseRepeat(bases, index, minRepeatCount);

            if(repeat == null)
            {
                repeat = findDualBaseRepeat(bases, index, minRepeatCount);
            }

            if(repeat == null)
            {
                repeat = findMultiBaseRepeat(bases, index, THREE_LENGTH, minRepeatCount);
            }

            if(repeat == null)
            {
                repeat = findMultiBaseRepeat(bases, index, FOUR_LENGTH, minRepeatCount);
            }

            if(repeat == null)
            {
                repeat = findMultiBaseRepeat(bases, index, FIVE_LENGTH, minRepeatCount);
            }

            if(repeat == null)
            {
                repeat = findDualDualRepeat(bases, index, minRepeatCount);
            }

            if(repeat != null)
            {
                if(repeats == null)
                    repeats = Lists.newArrayList(repeat);
                else
                    repeats.add(repeat);

                index += repeat.totalLength();
            }
            else
            {
                ++index;
            }
        }

        return repeats != null ? repeats : Collections.emptyList();
    }

    public static List<RepeatInfo> findRepeats(final byte[] bases, int startOffset, int endOffset)
    {
        if(startOffset == 0 && endOffset == 0)
            return findRepeats(bases);

        byte[] extensionBases = Arrays.copyOfRange(bases, startOffset, bases.length - endOffset);
        List<RepeatInfo> repeats = RepeatInfo.findRepeats(extensionBases);

        if(startOffset == 0)
            return repeats;

        // re-apply the offset skipped in the search
        return repeats.stream().map(x -> new RepeatInfo(x.Index + startOffset, x.Bases, x.Count)).collect(Collectors.toList());
    }

    // see Sage for a flexible repeat-length routine to find these
    private static final int DUAL_LENGTH = 2;
    private static final int THREE_LENGTH = 3;
    private static final int FOUR_LENGTH = 4;
    private static final int FIVE_LENGTH = 5;
    private static final int DUAL_DUAL_LENGTH = 4;

    public static RepeatInfo findSingleBaseRepeat(final byte[] bases, int index, final int minRepeatCount)
    {
        // the first base counts towards the repeat, so if length = 10, last index = 9, index cannot be higher than 9 - 4 + 1 = 6
        if(index + minRepeatCount - 1 >= bases.length)
            return null;

        int i = index + 1;
        while(i < bases.length)
        {
            if(bases[i] != bases[index])
                break;

            ++i;
        }

        int repeatCount = i - index;

        return repeatCount >= minRepeatCount ? new RepeatInfo(index, String.valueOf((char)bases[index]), repeatCount) : null;
    }

    public static RepeatInfo findDualBaseRepeat(final byte[] bases, int index, final int minRepeatCount)
    {
        if(index + minRepeatCount * 2 - 1 >= bases.length)
            return null;

        int i = index + DUAL_LENGTH;
        while(i < bases.length - 1)
        {
            if(bases[i] != bases[index] || bases[i + 1] != bases[index + 1])
                break;

            i += DUAL_LENGTH;
        }

        int repeatCount = (i - index) / DUAL_LENGTH;

        if(repeatCount < minRepeatCount)
            return null;

        String repeat = String.valueOf((char)bases[index]) + (char)bases[index + 1];
        return new RepeatInfo(index, repeat, repeatCount);
    }

    public static RepeatInfo findMultiBaseRepeat(final byte[] bases, int index, int repeatLength, final int minRepeatCount)
    {
        if(index + minRepeatCount * repeatLength - 1 >= bases.length)
            return null;

        int i = index + repeatLength;
        while(i < bases.length - (repeatLength - 1))
        {
            int matchedBases = 0;

            for(int j = 0; j < repeatLength; ++j)
            {
                if(bases[i + j] == bases[index + j])
                    ++matchedBases;
                else
                    break;
            }

            if(matchedBases != repeatLength)
                break;

            i += repeatLength;
        }

        int repeatCount = (i - index) / repeatLength;

        if(repeatCount < minRepeatCount)
            return null;

        String repeat = String.valueOf((char)bases[index]);

        for(int j = 1; j < repeatLength; ++j)
        {
            repeat += (char) bases[index + j];
        }

        return new RepeatInfo(index, repeat, repeatCount);
    }

    public static RepeatInfo findDualDualRepeat(final byte[] bases, int index, final int minRepeatCount)
    {
        if(index + minRepeatCount * DUAL_DUAL_LENGTH - 1 >= bases.length)
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

        if(repeatLength < minRepeatCount)
            return null;

        String repeat = String.valueOf((char)bases[index]) + (char)bases[index + 1] + (char)bases[index + 2] + (char)bases[index + 3];
        return new RepeatInfo(index, repeat, repeatLength);
    }

    public static int getRepeatCount(final Read read, final RepeatInfo repeatInfo, int readIndexStart, boolean searchForward)
    {
        return getRepeatCount(read.getBases(), repeatInfo.Bases, readIndexStart, searchForward);
    }

    public static int getRepeatCount(final byte[] bases, final String repeat, int readIndexStart, boolean searchForward)
    {
        // count how many instance of the repeat are in this read
        int repeatCount = 0;
        int repeatLength = repeat.length();
        int readIndex = readIndexStart;

        if(!searchForward)
            readIndex -= repeatLength - 1; // move to start of repeat

        if(readIndex < 0 || readIndex >= bases.length - repeatLength + 1)
            return -1;

        byte[] repeatBases = repeat.getBytes();

        while(readIndex >= 0 && readIndex < bases.length - repeatLength + 1)
        {
            for(int j = 0; j < repeatLength; ++j)
            {
                if(bases[readIndex + j] != repeatBases[j])
                    return repeatCount;
            }

            ++repeatCount;
            readIndex += searchForward ? repeatLength : -repeatLength;
        }

        return repeatCount;
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
                currentIndex += currentRepeat.totalLength();
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
                currentIndex += currentRepeat.totalLength();
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

    public static String repeatsAsString(final List<RepeatInfo> repeats)
    {
        if(repeats.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(" ");
        repeats.forEach(x -> sj.add(format("%dx%s", x.Count, x.Bases)));

        return format("%d %s", repeats.size(), sj);
    }
}
