package com.hartwig.hmftools.common.variant;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class Microhomology {
    private static final Logger LOGGER = LogManager.getLogger(Microhomology.class);

    private Microhomology() {
    }

    @NotNull
    public static String microhomologyAtDelete(int position, @NotNull final String sequence, @NotNull final String ref) {
        if (sequence.length() < position + ref.length()) {
            LOGGER.warn("Attempt to determine microhomology outside of sequence length");
            return Strings.EMPTY;
        }

        String result =
                commonPrefix(sequence.substring(position + 1), sequence.substring(position + 1 + ref.length() - 1), ref.length() - 1);

        for (int i = position; i >= 0; i--) {
            final String earlierPrefix = commonPrefix(sequence.substring(i), sequence.substring(i + ref.length() - 1), ref.length() - 1);
            if (earlierPrefix.length() > result.length()) {
                result = earlierPrefix;
            } else {
                return result;
            }
        }

        return result;
    }

    @NotNull
    public static String microhomologyAtInsert(int position, @NotNull final String sequence, @NotNull final String alt) {
        if (sequence.length() < position) {
            LOGGER.warn("Attempt to determine microhomology outside of sequence length");
            return Strings.EMPTY;
        }

        if (alt.contains(",")) {
            return Strings.EMPTY;
        }

        assert (sequence.charAt(position) == alt.charAt(0));
        final char ref = sequence.charAt(position);
        final String insert = alt.substring(1);
        final String preInsert = sequence.substring(0, position);
        final String postInsert = sequence.substring(position + 1);

        final String commonPrefix = commonPrefix(insert, postInsert);
        final String commonSuffixWithRef = commonSuffix(preInsert + ref, insert);

        return commonPrefix.length() > commonSuffixWithRef.length() ? commonPrefix : commonSuffixWithRef;
    }

    @VisibleForTesting
    @NotNull
    static String commonSuffix(@NotNull final String first, @NotNull final String second) {
        int minLength = Math.min(second.length(), first.length());
        if (minLength == 0) {
            return Strings.EMPTY;
        }

        for (int i = 0; i < minLength; i++) {
            if (first.charAt(first.length() - 1 - i) != second.charAt(second.length() - 1 - i)) {
                return i == 0 ? Strings.EMPTY : first.substring(first.length() - i);
            }
        }

        return first.substring(first.length() - minLength);
    }

    @NotNull
    private static String commonPrefix(@NotNull final String first, @NotNull final String second) {
        return commonPrefix(first, second, first.length());
    }

    @VisibleForTesting
    @NotNull
    static String commonPrefix(@NotNull final String first, @NotNull final String second, int maxLength) {
        int minLength = Math.min(maxLength, Math.min(second.length(), first.length()));
        if (minLength == 0) {
            return Strings.EMPTY;
        }

        for (int i = 0; i < minLength; i++) {
            if (first.charAt(i) != second.charAt(i)) {
                return i == 0 ? Strings.EMPTY : first.substring(0, i);
            }
        }
        return first.substring(0, minLength);
    }
}
