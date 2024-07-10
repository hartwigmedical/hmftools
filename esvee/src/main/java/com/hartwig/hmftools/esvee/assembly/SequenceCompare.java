package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

public final class SequenceCompare
{
    private static final int MISMATCH_RECOVERY_BASE_COUNT = 4;

    public static int calcReadSequenceMismatches(
            final boolean isForward, final byte[] extensionBases, final byte[] extensionQuals, final List<RepeatInfo> extensionRepeats,
            final Read read, final int readJunctionIndex, final int maxMismatches)
    {
        int readStartIndex = isForward ? readJunctionIndex : 0;
        int readEndIndex = isForward ? read.basesLength() - 1 : readJunctionIndex;

        // for -ve orientations, if extension sequence length = 10, with 0-8 being soft-clip and 9 being the first ref and junction index
        // and the read has 5 bases of soft-clip then read's start index will be 0 -> 4 + 1 = 5
        // so the comparison offset in the extension sequence is
        int extSeqReadStartIndex = isForward ? 0 : extensionBases.length - 1 - readJunctionIndex;

        byte[] readExtensionBases = subsetArray(read.getBases(), readStartIndex, readEndIndex);
        byte[] readExtensionQuals = subsetArray(read.getBaseQuality(), readStartIndex, readEndIndex);
        List<RepeatInfo> readRepeats = RepeatInfo.findRepeats(readExtensionBases);

        return compareSequences(
                extensionBases, extensionQuals, extSeqReadStartIndex, extensionBases.length - 1, extensionRepeats,
                readExtensionBases, readExtensionQuals, 0, readExtensionBases.length - 1,
                readRepeats != null ? readRepeats : Collections.emptyList(), maxMismatches);
    }

    public static boolean matchedAssemblySequences(final JunctionAssembly first, final JunctionAssembly second)
    {
        int firstIndexStart;
        int firstIndexEnd;
        int secondIndexStart;
        int secondIndexEnd;

        int junctionDiff = first.junction().Position - second.junction().Position;

        if(first.junction().isForward())
        {
            // where the junction position diff, take the index into the ref bases
            if(junctionDiff > 0)
            {
                firstIndexStart = first.junctionIndex() - junctionDiff;
                secondIndexStart = second.junctionIndex();
            }
            else
            {
                firstIndexStart = first.junctionIndex();
                secondIndexStart = second.junctionIndex() + junctionDiff;
            }

            int minDistanceFromJunction = min(first.upperDistanceFromJunction(), second.upperDistanceFromJunction());
            firstIndexEnd = firstIndexStart + minDistanceFromJunction;
            secondIndexEnd = secondIndexStart + minDistanceFromJunction;
        }
        else
        {
            int minDistanceFromJunction = min(first.lowerDistanceFromJunction(), second.lowerDistanceFromJunction());

            if(junctionDiff > 0)
            {
                firstIndexEnd = first.junctionIndex();
                secondIndexEnd = second.junctionIndex() + junctionDiff; // brings position back into the reference bases
            }
            else
            {
                firstIndexEnd = first.junctionIndex() - junctionDiff;
                secondIndexEnd = second.junctionIndex();
            }

            firstIndexStart = firstIndexEnd - minDistanceFromJunction;
            secondIndexStart = secondIndexEnd - minDistanceFromJunction;
        }

        if(firstIndexStart < 0 || secondIndexStart < 0)
            return false;

        firstIndexEnd = min(firstIndexEnd, first.bases().length - 1);
        secondIndexEnd = min(secondIndexEnd, second.bases().length - 1);

        if(firstIndexEnd <= firstIndexStart || secondIndexEnd <= secondIndexStart)
            return false;

        int mismatchCount = SequenceCompare.compareSequences(
                first.bases(), first.baseQuals(), firstIndexStart, firstIndexEnd, first.repeatInfo(),
                second.bases(), second.baseQuals(), secondIndexStart, secondIndexEnd, second.repeatInfo(), PRIMARY_ASSEMBLY_MERGE_MISMATCH);

        return mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH;
    }

    public static int compareSequences(
            final byte[] firstBases, final byte[] firstBaseQuals, int firstIndexStart, int firstIndexEnd, final List<RepeatInfo> firstRepeats,
            final byte[] secondBases, final byte[] secondBaseQuals, int secondIndexStart, int secondIndexEnd, final List<RepeatInfo> secondRepeats,
            int maxMismatches)
    {
        int mismatchCount = 0;

        int firstIndex = firstIndexStart;
        int secondIndex = secondIndexStart;

        int lastRepeatSkipBases = 0;

        while(firstIndex <= firstIndexEnd && secondIndex <= secondIndexEnd)
        {
            // the check for matching repeats is disabled since while logically a good idea, it is often throw out by prior mismatches
            // the solution may be to start the matching process from the other end
            if(basesMatch(
                    firstBases[firstIndex], secondBases[secondIndex], firstBaseQuals[firstIndex], secondBaseQuals[secondIndex],
                    LOW_BASE_QUAL_THRESHOLD))
            {
                ++firstIndex;
                ++secondIndex;
                continue;
            }

            // check if the mismatch is a single-base INDEL and if so skip it
            int recoverySkipBases = checkSkipShortIndel(firstBases, firstBaseQuals, firstIndex, secondBases, secondBaseQuals, secondIndex);

            if(recoverySkipBases != 0)
            {
                if(recoverySkipBases > 0)
                    firstIndex += recoverySkipBases;
                else
                    secondIndex += -recoverySkipBases;

                if(lastRepeatSkipBases == 0)
                {
                    if((recoverySkipBases == 1 && matchesRepeatStart(secondIndex, secondRepeats))
                    || (recoverySkipBases == -1 && matchesRepeatStart(firstIndex, firstRepeats)))
                    {
                        lastRepeatSkipBases = recoverySkipBases;
                    }
                }
                else if(lastRepeatSkipBases + recoverySkipBases == 0)
                {
                    lastRepeatSkipBases = 0;
                    continue;
                }

                ++mismatchCount;
                continue; // check the next base again
            }

            // check for a repeat diff - must be of the same type and just a different count
            int expectedRepeatBaseDiff = checkRepeatDifference(firstIndex, firstRepeats, secondIndex, secondRepeats);

            if(expectedRepeatBaseDiff != 0)
            {
                lastRepeatSkipBases = 0;

                // positive means the first's repeat was longer so move it ahead
                if(expectedRepeatBaseDiff > 0)
                    firstIndex += expectedRepeatBaseDiff;
                else
                    secondIndex += -expectedRepeatBaseDiff;

                ++mismatchCount;
                continue; // check the next base again
            }

            ++mismatchCount;

            if(maxMismatches >= 0 && mismatchCount > maxMismatches)
                return mismatchCount;

            ++firstIndex;
            ++secondIndex;
        }

        return mismatchCount;
    }

    private static int checkRepeatDifference(
            int firstIndex, final List<RepeatInfo> firstRepeats, int secondIndex, final List<RepeatInfo> secondRepeats)
    {
        // look for matching repeats, and if found return the expected difference in bases if the repeats have diff counts

        // example:
        // first index = 10, repeat length = 2, count = 5, so goes from 10-19
        // second index = 15, repeat length = 2, count = 7, so goes from 15-28
        // so would expect the repeat mismatch to occur at the base after either repeat ends

        RepeatInfo firstRepeat = findRepeat(firstIndex, firstRepeats);
        RepeatInfo secondRepeat = findRepeat(secondIndex, secondRepeats);

        if(firstRepeat == null || secondRepeat == null)
            return 0;

        if(!firstRepeat.matchesType(secondRepeat))
            return 0;

        if(firstRepeat.Count == secondRepeat.Count)
            return 0;

        return (firstRepeat.Count - secondRepeat.Count) * firstRepeat.Bases.length();
    }

    private static RepeatInfo findRepeat(int currentIndex, final List<RepeatInfo> repeats)
    {
        // find any repeat which covers or ends just prior to the current index
        for(RepeatInfo repeat: repeats)
        {
            if(repeat.postRepeatIndex() == currentIndex || (repeat.Index < currentIndex && repeat.postRepeatIndex() > currentIndex))
                return repeat;
        }

        return null;
    }

    private static boolean matchesRepeatStart(int currentIndex, final List<RepeatInfo> repeats)
    {
        // find any repeat which starts at the current index or next
        return repeats.stream().anyMatch(x -> x.Index == currentIndex || x.Index == currentIndex + 1);
    }

    private static int checkSkipShortIndel(
            final byte[] firstBases, final byte[] firstBaseQuals, final int firstIndex,
            final byte[] secondBases, final byte[] secondBaseQuals, final int secondIndex)
    {
        // test up to 2 skipped INDEL bases on each sequence, moving ahead 1 then 2 bases on the first sequence
        for(int diff = 1; diff <= 2; ++diff)
        {
            if(canRecoverMatch(
                    firstBases, firstBaseQuals, firstIndex + diff, secondBases, secondBaseQuals, secondIndex))
            {
                return diff;
            }
        }

        // then the same on the second sequence
        for(int diff = 1; diff <= 2; ++diff)
        {
            if(canRecoverMatch(
                    firstBases, firstBaseQuals, firstIndex, secondBases, secondBaseQuals, secondIndex + diff))
            {
                return -diff;
            }
        }

        return 0;
    }

    private static boolean canRecoverMatch(
            final byte[] firstBases, final byte[] firstBaseQuals, final int firstIndexStart,
            final byte[] secondBases, final byte[] secondBaseQuals, final int secondIndexStart)
    {
        int firstIndex = firstIndexStart;
        int secondIndex = secondIndexStart;
        int firstEndIndex = firstIndexStart + MISMATCH_RECOVERY_BASE_COUNT;
        int secondEndIndex = secondIndexStart + MISMATCH_RECOVERY_BASE_COUNT;

        if(firstEndIndex > firstBases.length || secondEndIndex > secondBases.length)
            return false;

        for(; firstIndex < firstEndIndex && secondIndex < secondEndIndex; ++firstIndex, ++secondIndex)
        {
            if(!basesMatch(
                    firstBases[firstIndex], secondBases[secondIndex], firstBaseQuals[firstIndex], secondBaseQuals[secondIndex],
                    LOW_BASE_QUAL_THRESHOLD))
            {
                return false;
            }
        }

        return true;
    }

    public static List<SequenceDiffInfo> getSequenceMismatchInfo(
            final byte[] firstBases, final byte[] firstBaseQuals, int firstIndexStart, int firstIndexEnd, final List<RepeatInfo> firstRepeats,
            final byte[] secondBases, final byte[] secondBaseQuals, int secondIndexStart, int secondIndexEnd, final List<RepeatInfo> secondRepeats)
    {
        List<SequenceDiffInfo> diffs = Lists.newArrayList();

        int firstIndex = firstIndexStart;
        int secondIndex = secondIndexStart;

        while(firstIndex <= firstIndexEnd && secondIndex <= secondIndexEnd)
        {
            // the check for matching repeats is disabled since while logically a good idea, it is often throw out by prior mismatches
            // the solution may be to start the matching process from the other end
            if(basesMatch(
                    firstBases[firstIndex], secondBases[secondIndex], firstBaseQuals[firstIndex], secondBaseQuals[secondIndex],
                    LOW_BASE_QUAL_THRESHOLD))
            {
                ++firstIndex;
                ++secondIndex;
                continue;
            }

            // check if the mismatch is a single-base INDEL and if so skip it
            int recoverySkipBases = checkSkipShortIndel(firstBases, firstBaseQuals, firstIndex, secondBases, secondBaseQuals, secondIndex);

            if(recoverySkipBases != 0)
            {
                if(recoverySkipBases > 0)
                    firstIndex += recoverySkipBases;
                else
                    secondIndex += -recoverySkipBases;

                boolean isFirst = recoverySkipBases > 0;

                diffs.add(new SequenceDiffInfo(
                        isFirst ? firstIndex : secondIndex, SequenceDiffType.INDEL, recoverySkipBases));

                continue; // check the next base again
            }

            // check for a repeat diff - must be of the same type and just a different count
            int expectedRepeatBaseDiff = checkRepeatDifference(firstIndex, firstRepeats, secondIndex, secondRepeats);

            if(expectedRepeatBaseDiff != 0)
            {
                // positive means the first's repeat was longer so move it ahead
                if(expectedRepeatBaseDiff > 0)
                    firstIndex += expectedRepeatBaseDiff;
                else
                    secondIndex += -expectedRepeatBaseDiff;

                boolean isFirst = expectedRepeatBaseDiff > 0;

                diffs.add(new SequenceDiffInfo(
                        isFirst ? firstIndex : secondIndex, SequenceDiffType.REPEAT, expectedRepeatBaseDiff));

                continue; // check the next base again
            }

            diffs.add(new SequenceDiffInfo(firstIndex, SequenceDiffType.BASE, 1));

            ++firstIndex;
            ++secondIndex;
        }

        return diffs;
    }
}
