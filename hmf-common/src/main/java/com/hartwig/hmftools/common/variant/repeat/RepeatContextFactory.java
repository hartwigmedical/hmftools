package com.hartwig.hmftools.common.variant.repeat;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class RepeatContextFactory {

    private static final int MIN_COUNT = 2;
    private static final int MAX_LENGTH = 10;

    private static final Logger LOGGER = LogManager.getLogger(RepeatContextFactory.class);

    private RepeatContextFactory() {
    }

    @NotNull
    public static Optional<RepeatContext> repeats(int index, final byte[] readSequence) {

        final List<RepeatContext> repeatContexts = Lists.newArrayList();

        if (readSequence.length >= index) {

            for (int repeatStartIndex = Math.max(0, index - MAX_LENGTH); repeatStartIndex <= index; repeatStartIndex++) {
                for (int repeatEndIndex = index;
                        repeatEndIndex <= Math.min(readSequence.length, repeatStartIndex + MAX_LENGTH); repeatEndIndex++) {

                    int repeatLength = repeatEndIndex - repeatStartIndex + 1;
                    int forwardsCount = forwardRepeats(repeatStartIndex, repeatLength, readSequence);
                    int backwardsCount = backwardRepeats(repeatStartIndex, repeatLength, readSequence);
                    if (forwardsCount + backwardsCount > MIN_COUNT) {
                        repeatContexts.add(new RepeatContext(readSequence, repeatStartIndex, repeatLength, forwardsCount, backwardsCount));
                    }
                }
            }
        } else {
            LOGGER.warn("Repeats requested outside of sequence length");
        }

        repeatContexts.sort(Comparator.comparingInt(RepeatContext::count).reversed());
        return repeatContexts.isEmpty() ? Optional.empty() : Optional.of(repeatContexts.get(0));
    }

    @NotNull
    public static Optional<RepeatContext> repeats(int index, @NotNull final String sequence) {
        return repeats(index, sequence.getBytes());
    }

    @VisibleForTesting
    static int forwardRepeats(int index, int repeatLength, final byte[] readSequence) {
        for (int count = 1; ; count++) {
            if (!match(index, repeatLength, index + count * repeatLength, readSequence)) {
                return count;
            }
        }
    }

    @VisibleForTesting
    static int backwardRepeats(int index, int repeatLength, final byte[] readSequence) {
        for (int count = 1; ; count++) {
            if (!match(index, repeatLength, index - count * repeatLength, readSequence)) {
                return count - 1;
            }
        }
    }

    @VisibleForTesting
    static boolean match(int repeatIndex, int repeatLength, int readIndex, byte[] readSequence) {
        if (readIndex + repeatLength > readSequence.length || readIndex < 0 || repeatIndex + repeatLength > readSequence.length) {
            return false;
        }

        for (int i = 0; i < repeatLength; i++) {
            if (readSequence[repeatIndex + i] != readSequence[readIndex + i]) {
                return false;
            }
        }

        return true;
    }

}
