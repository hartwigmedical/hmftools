package com.hartwig.hmftools.sage.old;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;

public final class RepeatContextFactory
{
    private static final int MIN_COUNT = 2;
    private static final int MAX_LENGTH = 10;

    public static RepeatContext findRepeat(int index, final byte[] readBases, final int minCount, final int maxLength)
    {
        if(index > readBases.length)
            return null;

        RepeatContext maxRepeatCount = null;

        for(int repeatStartIndex = max(0, index - maxLength); repeatStartIndex <= index; repeatStartIndex++)
        {
            for(int repeatEndIndex = index; repeatEndIndex <= min(readBases.length, repeatStartIndex + maxLength); repeatEndIndex++)
            {
                int repeatLength = repeatEndIndex - repeatStartIndex + 1;
                int forwardsCount = forwardRepeats(repeatStartIndex, repeatLength, readBases);
                int backwardsCount = backwardRepeats(repeatStartIndex, repeatLength, readBases);

                if(forwardsCount + backwardsCount >= minCount)
                {
                    int startIndex = repeatStartIndex - backwardsCount * repeatLength;
                    int endIndex = repeatStartIndex + forwardsCount * repeatLength - 1;
                    int additionalBasesAtEnd = matchingBasesFromLeft(repeatStartIndex, repeatLength, endIndex + 1, readBases);

                    if(maxRepeatCount == null || forwardsCount + backwardsCount > maxRepeatCount.count())
                    {
                        maxRepeatCount = new RepeatContext(
                                readBases, repeatStartIndex, startIndex, endIndex + additionalBasesAtEnd,
                                repeatLength, forwardsCount, backwardsCount);
                    }
                }
            }
        }

        return maxRepeatCount;
    }

    public static Optional<RepeatContext> repeats(int index, final byte[] readSequence)
    {
        RepeatContext repeatContext = findRepeat(index, readSequence, MIN_COUNT, MAX_LENGTH);
        return repeatContext != null ? Optional.of(repeatContext) : Optional.empty();
    }

    public static Optional<RepeatContext> repeats(int index, final String sequence)
    {
        return repeats(index, sequence.getBytes());
    }

    @VisibleForTesting
    public static int forwardRepeats(int index, int repeatLength, final byte[] readSequence)
    {
        for(int count = 1; ; count++)
        {
            if(!match(index, repeatLength, index + count * repeatLength, readSequence))
            {
                return count;
            }
        }
    }

    public static int backwardRepeats(int index, int repeatLength, final byte[] readSequence)
    {
        for(int count = 1; ; count++)
        {
            if(!match(index, repeatLength, index - count * repeatLength, readSequence))
            {
                return count - 1;
            }
        }
    }

    @VisibleForTesting
    public static boolean match(int repeatIndex, int repeatLength, int readIndex, byte[] readSequence)
    {
        return matchingBasesFromLeft(repeatIndex, repeatLength, readIndex, readSequence) == repeatLength;
    }

    private static int matchingBasesFromLeft(int repeatIndex, int repeatLength, int readIndex, byte[] readSequence)
    {
        for(int i = 0; i < repeatLength; i++)
        {
            if(outOfBounds(repeatIndex + i, readSequence) || outOfBounds(readIndex + i, readSequence)
                    || readSequence[repeatIndex + i] != readSequence[readIndex + i])
            {
                return i;
            }
        }

        return repeatLength;
    }

    private static boolean outOfBounds(int index, byte[] sequence)
    {
        return index < 0 || index >= sequence.length;
    }
}
