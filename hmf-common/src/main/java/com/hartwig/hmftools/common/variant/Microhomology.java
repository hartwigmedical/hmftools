package com.hartwig.hmftools.common.variant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class Microhomology {

    private Microhomology() {
    }

    @NotNull
    static String microhomology(int position, @NotNull final String sequence, @NotNull final String ref) {
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
    private static String commonPrefix(String normal, String tumor, int maxLength) {
        int minLength = Math.min(maxLength, Math.min(tumor.length(), normal.length()));
        if (minLength == 0) {
            return Strings.EMPTY;
        }

        for (int i = 0; i < minLength; i++) {
            if (normal.charAt(i) != tumor.charAt(i)) {
                return i == 0 ? Strings.EMPTY : normal.substring(0, i);
            }
        }
        return normal.substring(0, minLength);
    }
}
