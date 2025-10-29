package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REPEAT_2_DIFF_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REPEAT_3_DIFF_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
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

        // check for a simple match of one assembly's extension bases present in the other and straddling the other's junction
        String firstFullSequence = first.formFullSequence();
        String firstOverlapSequence = firstFullSequence.substring(firstIndexStart, firstIndexEnd + 1);

        String secondFullSequence = second.formFullSequence();

        int secondExtBaseIndex = secondFullSequence.indexOf(firstOverlapSequence);

        if(secondExtBaseIndex >= 0)
        {
            int secondExtBaseEndIndex = secondExtBaseIndex + firstOverlapSequence.length() - 1;
            if(secondExtBaseIndex < second.junctionIndex() && secondExtBaseEndIndex > second.junctionIndex())
                return true;
        }

        String secondOverlapSequence = secondFullSequence.substring(secondIndexStart, secondIndexEnd + 1);

        int firstExtBaseIndex = firstFullSequence.indexOf(secondOverlapSequence);

        if(firstExtBaseIndex >= 0)
        {
            int firstExtBaseEndIndex = firstExtBaseIndex + secondOverlapSequence.length() - 1;
            if(firstExtBaseIndex < first.junctionIndex() && firstExtBaseEndIndex > first.junctionIndex())
                return true;
        }

        int mismatchCount = compareSequences(
                first.bases(), first.baseQuals(), firstIndexStart, firstIndexEnd, first.repeatInfo(),
                second.bases(), second.baseQuals(), secondIndexStart, secondIndexEnd, second.repeatInfo(), PRIMARY_ASSEMBLY_MERGE_MISMATCH);

        if(mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH)
            return true;

        /* NOTE: reconsider this once the compareSequences has improved repeat handling logic

        // try subsequence matches to anchor the full comparison
        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;

        int firstSubSeqIndex = 0;

        int subSeqIterations = (int)floor(firstFullSequence.length() / subsequenceLength);

        for(int i = 0; i < subSeqIterations; ++i)
        {
            firstSubSeqIndex = i * subsequenceLength;
            int matchSeqEndIndex = firstSubSeqIndex + subsequenceLength;

            if(matchSeqEndIndex >= firstFullSequence.length())
                break;

            String assemblySubSequence = firstFullSequence.substring(firstSubSeqIndex, firstSubSeqIndex + subsequenceLength);

            int secondSubSeqIndex = secondFullSequence.indexOf(assemblySubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            int minLowerOverlap = min(firstSubSeqIndex, secondSubSeqIndex);
            int minUpperOverlap = min(firstFullSequence.length() - firstSubSeqIndex - 1, secondFullSequence.length() - secondSubSeqIndex - 1);

            if(minLowerOverlap + minUpperOverlap < PROXIMATE_JUNCTION_DISTANCE)
                continue;

            firstIndexStart = firstSubSeqIndex - minLowerOverlap;
            secondIndexStart = secondSubSeqIndex - minLowerOverlap;

            firstIndexEnd = min(firstSubSeqIndex + minUpperOverlap, first.bases().length - 1);
            secondIndexEnd = min(secondSubSeqIndex + minUpperOverlap, second.bases().length - 1);

            mismatchCount = compareSequences(
                    first.bases(), first.baseQuals(), firstIndexStart, firstIndexEnd, first.repeatInfo(),
                    second.bases(), second.baseQuals(), secondIndexStart, secondIndexEnd, second.repeatInfo(), PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH)
                return true;
        }
        */

        return false;
    }

    public static int compareSequences(
            final byte[] firstBases, final byte[] firstBaseQuals, int firstIndexStart, int firstIndexEnd, final List<RepeatInfo> firstRepeats,
            final byte[] secondBases, final byte[] secondBaseQuals, int secondIndexStart, int secondIndexEnd, final List<RepeatInfo> secondRepeats,
            int maxMismatches)
    {
        return compareSequencesMismatchSpec(
                firstBases, firstBaseQuals, firstIndexStart, firstIndexEnd, firstRepeats, secondBases, secondBaseQuals, secondIndexStart,
                secondIndexEnd, secondRepeats, maxMismatches, true, true);
    }

    private static int compareSequencesMismatchSpec(
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
        if(repeatCount >= REPEAT_3_DIFF_COUNT)
            return 3;

        if(repeatCount >= REPEAT_2_DIFF_COUNT)
            return 2;

        return 1;
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
}
