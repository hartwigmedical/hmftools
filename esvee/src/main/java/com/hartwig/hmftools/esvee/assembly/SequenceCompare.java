package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REPEAT_2_DIFF_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

public final class SequenceCompare
{
    private static final int MISMATCH_RECOVERY_BASE_COUNT = 3;

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

        int mismatchCount = compareSequences(
                first.bases(), first.baseQuals(), firstIndexStart, firstIndexEnd, first.repeatInfo(),
                second.bases(), second.baseQuals(), secondIndexStart, secondIndexEnd, second.repeatInfo(), PRIMARY_ASSEMBLY_MERGE_MISMATCH);

        return mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH;
    }

    public static int compareSequences(
            final byte[] firstBases, final byte[] firstBaseQuals, int firstIndexStart, int firstIndexEnd, final List<RepeatInfo> firstRepeats,
            final byte[] secondBases, final byte[] secondBaseQuals, int secondIndexStart, int secondIndexEnd, final List<RepeatInfo> secondRepeats,
            int maxMismatches)
    {
        return compareSequences(
                firstBases, firstBaseQuals, firstIndexStart, firstIndexEnd, firstRepeats, secondBases, secondBaseQuals, secondIndexStart,
                secondIndexEnd, secondRepeats, maxMismatches, true, true);
    }

    public static int compareSequences(
            final byte[] firstBases, final byte[] firstBaseQuals, int firstIndexStart, int firstIndexEnd, final List<RepeatInfo> firstRepeats,
            final byte[] secondBases, final byte[] secondBaseQuals, int secondIndexStart, int secondIndexEnd, final List<RepeatInfo> secondRepeats,
            int maxMismatches, boolean checkForwards, boolean repeatsAreMismatches)
    {
        if(firstIndexStart < 0 || firstIndexEnd >= firstBases.length || firstIndexEnd >= firstBaseQuals.length
        || secondIndexStart < 0 || secondIndexEnd >= secondBases.length || secondIndexEnd >= secondBaseQuals.length)
        {
            return maxMismatches < 0 ? maxMismatches : maxMismatches + 1; // invalid
        }

        int mismatchCount = 0;

        int firstIndex = checkForwards ? firstIndexStart : firstIndexEnd;
        int secondIndex = checkForwards ? secondIndexStart : secondIndexEnd;

        int lastRepeatSkipBases = 0;

        while(firstIndex >= 0 && firstIndex <= firstIndexEnd && secondIndex >= 0 && secondIndex <= secondIndexEnd)
        {
            // the check for matching repeats is disabled since while logically a good idea, it is often throw out by prior mismatches
            // the solution may be to start the matching process from the other end
            if(basesMatch(firstBases[firstIndex], secondBases[secondIndex], firstBaseQuals[firstIndex], secondBaseQuals[secondIndex]))
            {
                firstIndex += checkForwards ? 1 : -1;
                secondIndex += checkForwards ? 1 : -1;
                continue;
            }

            // check if the mismatch is a single-base INDEL and if so skip it
            int recoverySkipBases = checkSkipShortIndel(
                    firstBases, firstBaseQuals, firstIndex, secondBases, secondBaseQuals, secondIndex, checkForwards);

            if(recoverySkipBases != 0)
            {
                if(recoverySkipBases > 0)
                    firstIndex += checkForwards ? recoverySkipBases : -recoverySkipBases;
                else
                    secondIndex += checkForwards ? -recoverySkipBases : recoverySkipBases;

                if(lastRepeatSkipBases == 0)
                {
                    if((recoverySkipBases == 1 && matchesRepeat(secondIndex, secondRepeats, checkForwards))
                    || (recoverySkipBases == -1 && matchesRepeat(firstIndex, firstRepeats, checkForwards)))
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
            int expectedRepeatBaseDiff = checkRepeatDifference(firstIndex, firstRepeats, secondIndex, secondRepeats, checkForwards);

            if(expectedRepeatBaseDiff != 0)
            {
                lastRepeatSkipBases = 0;

                // positive means the first's repeat was longer so move it ahead
                if(expectedRepeatBaseDiff > 0)
                    firstIndex += checkForwards ? expectedRepeatBaseDiff : expectedRepeatBaseDiff;
                else
                    secondIndex += checkForwards ? -expectedRepeatBaseDiff : expectedRepeatBaseDiff;

                if(repeatsAreMismatches)
                    ++mismatchCount;

                continue; // check the next base again
            }

            ++mismatchCount;

            if(maxMismatches >= 0 && mismatchCount > maxMismatches)
                return mismatchCount;

            firstIndex += checkForwards ? 1 : -1;
            secondIndex += checkForwards ? 1 : -1;
        }

        return mismatchCount;
    }

    public static int permittedRepeatCount(final int repeatCount)
    {
        return repeatCount < REPEAT_2_DIFF_COUNT ? 1 : 2;
    }

    private static int checkRepeatDifference(
            int firstIndex, final List<RepeatInfo> firstRepeats, int secondIndex, final List<RepeatInfo> secondRepeats, boolean checkForwards)
    {
        // look for matching repeats, and if found return the expected difference in bases if the repeats have diff counts

        // example:
        // first index = 10, repeat length = 2, count = 5, so goes from 10-19
        // second index = 15, repeat length = 2, count = 7, so goes from 15-28
        // so would expect the repeat mismatch to occur at the base after either repeat ends

        RepeatInfo firstRepeat = findRepeat(firstIndex, firstRepeats, checkForwards);
        RepeatInfo secondRepeat = findRepeat(secondIndex, secondRepeats, checkForwards);

        if(firstRepeat == null && secondRepeat == null)
            return 0;

        if(firstRepeat != null && secondRepeat != null)
        {
            if(!firstRepeat.matchesType(secondRepeat))
            {
                // consider the longer repeat if it contains the shorter and has the minimum count
                // how should this make use of max repeat count differences?
                if(firstRepeat.Bases.contains(secondRepeat.Bases) && firstRepeat.Count == 2)
                {
                    int ratio = firstRepeat.baseLength() / secondRepeat.baseLength();

                    if(secondRepeat.Count != firstRepeat.Count * ratio)
                        return 0;

                    return firstRepeat.Bases.length();
                }
                else if(secondRepeat.Bases.contains(firstRepeat.Bases) && secondRepeat.Count == 2)
                {
                    int ratio = secondRepeat.baseLength() / firstRepeat.baseLength();

                    if(firstRepeat.Count != secondRepeat.Count * ratio)
                        return 0;

                    return secondRepeat.Bases.length();
                }
                return 0;
            }

            if(firstRepeat.Count == secondRepeat.Count)
                return 0;

            int permittedRepeatDiff = permittedRepeatCount(max(firstRepeat.Count, secondRepeat.Count));

            if(abs(firstRepeat.Count - secondRepeat.Count) > permittedRepeatDiff)
                return 0;

            return (firstRepeat.Count - secondRepeat.Count) * firstRepeat.Bases.length();
        }

        /*
        // also permitted is a repeat in only one of the sequences at the limit of min repeat count
        // NOTE: this would require checking the other sequence's bases at this location to ensure they maatched the repeat
        if(firstRepeat != null && secondRepeat == null)
        {
            if(firstRepeat.Count == minRequiredRepeatCount(firstRepeat.baseLength()))
                return firstRepeat.baseLength();
        }
        else if(firstRepeat == null && secondRepeat != null)
        {
            if(secondRepeat.Count == minRequiredRepeatCount(secondRepeat.baseLength()))
                return -secondRepeat.baseLength();
        }
        */

        return 0; // no valid repeat adjustment
    }

    private static RepeatInfo findRepeat(int currentIndex, final List<RepeatInfo> repeats, boolean checkForwards)
    {
        // find any repeat which covers or ends just prior to the current index
        for(RepeatInfo repeat: repeats)
        {
            if(repeat.Index < currentIndex && repeat.postRepeatIndex() > currentIndex)
                return repeat;

            if(checkForwards)
            {
                if(repeat.postRepeatIndex() == currentIndex)
                    return repeat;
            }
            else
            {
                if(repeat.Index == currentIndex)
                    return repeat;
            }
        }

        return null;
    }

    private static boolean matchesRepeat(int currentIndex, final List<RepeatInfo> repeats, boolean checkForwards)
    {
        return checkForwards ? matchesRepeatStart(currentIndex, repeats) : matchesRepeatEnd(currentIndex, repeats);
    }

    private static boolean matchesRepeatStart(int currentIndex, final List<RepeatInfo> repeats)
    {
        // find any repeat which starts at the current index or next
        return repeats.stream().anyMatch(x -> x.Index == currentIndex || x.Index == currentIndex + 1);
    }

    private static boolean matchesRepeatEnd(int currentIndex, final List<RepeatInfo> repeats)
    {
        // find any repeat which starts at the current index or next
        return repeats.stream().anyMatch(x -> x.postRepeatIndex() == currentIndex || x.postRepeatIndex() == currentIndex - 1);
    }

    private static int checkSkipShortIndel(
            final byte[] firstBases, final byte[] firstBaseQuals, final int firstIndex,
            final byte[] secondBases, final byte[] secondBaseQuals, final int secondIndex, boolean checkForwards)
    {
        // test up to 2 skipped INDEL bases on each sequence, moving ahead 1 then 2 bases on the first sequence
        for(int diff = 1; diff <= 2; ++diff)
        {
            int adjustedFirstIndex = firstIndex + (checkForwards ? diff : -diff);

            if(canRecoverMatch(
                    firstBases, firstBaseQuals, adjustedFirstIndex, secondBases, secondBaseQuals, secondIndex, checkForwards))
            {
                return diff;
            }
        }

        // then the same on the second sequence
        for(int diff = 1; diff <= 2; ++diff)
        {
            int adjustedSecondIndex = secondIndex + (checkForwards ? diff : -diff);

            if(canRecoverMatch(
                    firstBases, firstBaseQuals, firstIndex, secondBases, secondBaseQuals, adjustedSecondIndex, checkForwards))
            {
                return -diff;
            }
        }

        return 0;
    }

    private static boolean canRecoverMatch(
            final byte[] firstBases, final byte[] firstBaseQuals, final int firstCurrentIndex,
            final byte[] secondBases, final byte[] secondBaseQuals, final int secondCurrentStart, boolean checkForwards)
    {
        int firstIndex = firstCurrentIndex;
        int secondIndex = secondCurrentStart;

        int firstIndexEnd, secondIndexEnd;

        if(checkForwards)
        {
            firstIndexEnd = firstCurrentIndex + MISMATCH_RECOVERY_BASE_COUNT;
            secondIndexEnd = secondCurrentStart + MISMATCH_RECOVERY_BASE_COUNT;

            if(firstIndexEnd >= firstBases.length || secondIndexEnd >= secondBases.length)
                return false;
        }
        else
        {
            firstIndexEnd = firstCurrentIndex - MISMATCH_RECOVERY_BASE_COUNT;
            secondIndexEnd = secondCurrentStart - MISMATCH_RECOVERY_BASE_COUNT;

            if(firstIndexEnd < 0 || secondIndexEnd < 0)
                return false;
        }

        while(true)
        {
            if(!basesMatch(firstBases[firstIndex], secondBases[secondIndex], firstBaseQuals[firstIndex], secondBaseQuals[secondIndex]))
            {
                return false;
            }

            if(firstIndex == firstIndexEnd && secondIndex == secondIndexEnd)
                break;

            firstIndex += checkForwards ? 1 : -1;
            secondIndex += checkForwards ? 1 : -1;
        }

        return true;
    }

    // unused matching routines
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

        if(extSeqReadStartIndex < 0)
        {
            readStartIndex += -extSeqReadStartIndex;
            extSeqReadStartIndex = 0;
        }

        byte[] readExtensionBases = subsetArray(read.getBases(), readStartIndex, readEndIndex);
        byte[] readExtensionQuals = subsetArray(read.getBaseQuality(), readStartIndex, readEndIndex);
        List<RepeatInfo> readRepeats = RepeatInfo.findRepeats(readExtensionBases);

        // run the comparison heading out from the junction
        return compareSequences(
                extensionBases, extensionQuals, extSeqReadStartIndex, extensionBases.length - 1, extensionRepeats,
                readExtensionBases, readExtensionQuals, 0, readExtensionBases.length - 1,
                readRepeats != null ? readRepeats : Collections.emptyList(), maxMismatches, isForward, false);
    }

    private static boolean hasCompatibleRepeats(
            int firstIndexStart, int firstIndexEnd, final List<RepeatInfo> firstRepeats,
            int secondIndexStart, int secondIndexEnd, final List<RepeatInfo> secondRepeats)
    {
        List<RepeatInfo> firstRepeatsInRange = firstRepeats.stream()
                .filter(x -> x.Index >= firstIndexStart && x.postRepeatIndex() <= firstIndexEnd).collect(Collectors.toList());

        List<RepeatInfo> secondRepeatsInRange = secondRepeats.stream()
                .filter(x -> x.Index >= secondIndexStart && x.postRepeatIndex() <= secondIndexEnd).collect(Collectors.toList());

        int minRepeatCount = min(firstRepeatsInRange.size(), secondRepeatsInRange.size());

        int firstIndex = 0;
        int matchedCount = 0;

        while(firstIndex < firstRepeatsInRange.size())
        {
            RepeatInfo firstRepeat = firstRepeatsInRange.get(firstIndex);

            RepeatInfo matchedRepeat = secondRepeatsInRange.stream()
                    .filter(x -> x.Bases.equals(firstRepeat.Bases) && x.Count == firstRepeat.Count).findFirst().orElse(null);

            if(matchedRepeat == null)
            {
                matchedRepeat = secondRepeatsInRange.stream()
                        .filter(x -> x.Bases.equals(firstRepeat.Bases) && abs(x.Count - firstRepeat.Count) <= 2).findFirst().orElse(null);
            }

            if(matchedRepeat != null)
            {
                ++matchedCount;
                firstRepeatsInRange.remove(firstIndex);
                secondRepeatsInRange.remove(matchedRepeat);
            }
            else
            {
                ++firstIndex;
            }
        }

        return minRepeatCount - matchedCount <= 2;
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
            if(basesMatch(firstBases[firstIndex], secondBases[secondIndex], firstBaseQuals[firstIndex], secondBaseQuals[secondIndex]))
            {
                ++firstIndex;
                ++secondIndex;
                continue;
            }

            // check if the mismatch is a single-base INDEL and if so skip it
            int recoverySkipBases = checkSkipShortIndel(firstBases, firstBaseQuals, firstIndex, secondBases, secondBaseQuals, secondIndex, true);

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
            int expectedRepeatBaseDiff = checkRepeatDifference(firstIndex, firstRepeats, secondIndex, secondRepeats, true);

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
