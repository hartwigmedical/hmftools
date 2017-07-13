package com.hartwig.hmftools.common.purple.repeat;

import java.util.Comparator;
import java.util.Map;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class RepeatContextFactory {

    private static final int MIN_COUNT = 4;
    private static final int MAX_LENGTH = 6;

    public static Optional<RepeatContext> repeats(int position, @NotNull final String refGenome, @NotNull final String ref, @NotNull final String alt) {
        if (isDelete(ref, alt)) {
            return repeats(position + 1, refGenome);
        }
        return Optional.empty();
    }

    private static boolean isDelete(@NotNull final String ref, @NotNull final String alt) {
        return alt.length() == 1 && ref.length() > 1;
    }

    private static boolean isSimpleIndel(@NotNull final String ref, @NotNull final String alt) {
        return ref.length() != alt.length() && (ref.length() == 1 || alt.length() == 1);
    }

    public static Optional<RepeatContext> repeats(int index, @NotNull final String sequence) {
        final Map<String, Integer> result = Maps.newHashMap();

        for (int start = Math.max(0, index - MAX_LENGTH); start <= index; start++) {

            final String prior = sequence.substring(0, start);
            final String post = sequence.substring(start);

            for (int end = index; end <= Math.min(sequence.length(), start + MAX_LENGTH); end++) {
                if (end != index) {

                    int count = 0;
                    final String bases = sequence.substring(Math.min(start, end), Math.max(start, end));

                    count += backwardRepeats(bases, prior);
                    count += forwardRepeats(bases, post);

                    if (count >= MIN_COUNT) {
                        result.merge(bases, count, Math::max);
                    }
                }
            }
        }
        return result.entrySet().stream().max(Comparator.comparingInt(Map.Entry::getValue)).map(RepeatContextFactory::create);
    }

    private static RepeatContext create(Map.Entry<String, Integer> entry) {
        return ImmutableRepeatContext.builder().sequence(entry.getKey()).count(entry.getValue()).build();
    }

    @VisibleForTesting
    static int forwardRepeats(@NotNull final String bases, @NotNull final String sequence) {
        int count = 0;
        int basesLength = bases.length();
        for (int j = 0; j < sequence.length() - basesLength + 1; j += basesLength) {
            final String subSequence = sequence.substring(j, j + basesLength);
            if (!subSequence.equals(bases)) {
                break;
            }
            count++;
        }

        return count;
    }

    @VisibleForTesting
    static int backwardRepeats(@NotNull final String bases, @NotNull final String sequence) {
        int count = 0;
        int basesLength = bases.length();

        for (int j = sequence.length() - basesLength; j >= 0; j -= basesLength) {
            final String subSequence = sequence.substring(j, j + basesLength);
            if (!subSequence.equals(bases)) {
                break;
            }
            count++;
        }

        return count;
    }

}
