package com.hartwig.hmftools.common.variant.repeat;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RepeatContextFactory_
{
    private static final int MIN_COUNT = 2;
    private static final int MAX_LENGTH = 10;

    private static final Logger LOGGER = LogManager.getLogger(RepeatContextFactory_.class);

    /**
     * Searches for a repeat ending starting between [index - MAX_LENGTH, index], and ending from [index, repeatStartIndex + MAX_LENGTH]
     * clipped via readSequence boundaries.
     *
     * If multiple repeates are found then return the one with the largest number of repeats. So chooses the "core" of the repeat.
     */
    @NotNull
    public static Optional<RepeatContext> repeats(int index, final byte[] readSequence)
    {
        if (index > readSequence.length)
        {
            LOGGER.warn("Repeats requested outside of sequence length");
            return Optional.empty();
        }

        List<RepeatContext> repeatContexts = Lists.newArrayList();
        for (int repeatStartIndex = max(0, index - MAX_LENGTH); repeatStartIndex <= index; ++repeatStartIndex)
        {
            for (int repeatEndIndex = index; repeatEndIndex <= min(readSequence.length, repeatStartIndex + MAX_LENGTH); ++repeatEndIndex)
            {
                Optional<RepeatContext> repeatContext = getRepeatContext(repeatStartIndex, repeatEndIndex, readSequence);
                if (repeatContext.isPresent())
                {
                    repeatContexts.add(repeatContext.get());
                }
            }
        }

        repeatContexts.sort(Comparator.comparingInt(RepeatContext::count).reversed());
        return repeatContexts.isEmpty() ? Optional.empty() : Optional.of(repeatContexts.get(0));
    }

    /**
     * Returns the ReadContext based on a repeat from repeatStartIndex to repeatEndIndex.
     * <p>
     * Has to have at least MIN_COUNT repeats.
     */
    @NotNull
    private static Optional<RepeatContext> getRepeatContext(int repeatStartIndex, int repeatEndIndex, final byte[] readSequence)
    {
        int repeatLength = repeatEndIndex - repeatStartIndex + 1;
        int forwardCount = forwardRepeats(repeatStartIndex, repeatLength, readSequence);
        int backwardCount = backwardRepeats(repeatStartIndex, repeatLength, readSequence);
        if (forwardCount + backwardCount < MIN_COUNT)
            return Optional.empty();

        int startIndex = repeatStartIndex - backwardCount * repeatLength;
        int endIndex = repeatStartIndex + forwardCount * repeatLength - 1;
        int extraBasesAtEnd = matchingBasesFromLeft(repeatStartIndex, repeatLength, endIndex + 1, readSequence);

        return Optional.of(new RepeatContext(
                readSequence,
                repeatStartIndex,
                startIndex,
                endIndex + extraBasesAtEnd,
                repeatLength,
                forwardCount,
                backwardCount
        ));
    }

    @NotNull
    public static Optional<RepeatContext> repeats(int index, @NotNull final String sequence)
    {
        return repeats(index, sequence.getBytes());
    }

    /**
     * Returns how many times the repeat starting at index repeats forwards.
     */
    @VisibleForTesting
    public static int forwardRepeats(int index, int repeatLength, final byte[] readSequence)
    {
        for (int repeatCount = 1; ; ++repeatCount)
        {
            if (!match(index, repeatLength, index + repeatCount * repeatLength, readSequence))
                return repeatCount;
        }
    }

    /**
     * Returns how many times the repeat starting at index repeats backwards (does not include initial repeat).
     */
    public static int backwardRepeats(int index, int repeatLength, final byte[] readSequence)
    {
        for (int repeatCount = 1; ; ++repeatCount)
        {
            if (!match(index, repeatLength, index - repeatCount * repeatLength, readSequence))
                return repeatCount - 1;
        }
    }

    /**
     * Returns whether repeatLength bases in readSequence starting from readIndex matches repeatLength bases starting from repeatIndex.
     */
    @VisibleForTesting
    static boolean match(int repeatIndex, int repeatLength, int readIndex, byte[] readSequence)
    {
        return matchingBasesFromLeft(repeatIndex, repeatLength, readIndex, readSequence) == repeatLength;
    }

    /**
     * Counts the number of bases in readSequence starting from readIndex matches the bases starting from repeatIndex up to repeatLength.
     */
    private static int matchingBasesFromLeft(int repeatIndex, int repeatLength, int readIndex, byte[] readSequence)
    {
        for (int i = 0; i < repeatLength; ++i)
        {
            if (outOfBounds(repeatIndex + i, readSequence) || outOfBounds(readIndex + i, readSequence) || readSequence[repeatIndex + i] != readSequence[readIndex + i])
                return i;
        }

        return repeatLength;
    }

    private static boolean outOfBounds(int index, byte[] sequence)
    {
        return index < 0 || index >= sequence.length;
    }
}
