package com.hartwig.hmftools.vicc.annotation;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;
import com.hartwig.hmftools.common.serve.classification.ExclusiveEventMatcher;

import org.jetbrains.annotations.NotNull;

class GeneRangeCodonClassifier implements EventMatcher {

    @NotNull
    public static EventMatcher create(@NotNull List<EventMatcher> excludingEventMatchers) {
        return new ExclusiveEventMatcher(excludingEventMatchers, new GeneRangeCodonClassifier());
    }

    private GeneRangeCodonClassifier() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        String proteinAnnotation = HotspotClassifier.extractProteinAnnotation(event);

        return isValidSingleCodonRange(proteinAnnotation);
    }

    private static boolean isValidSingleCodonRange(@NotNull String event) {
        // Feature codon ranges occasionally come with parentheses
        String strippedEvent = event.replace("(", "").replace(")", "");

        // Features are expected to look something like V600 (1 char - N digits)
        if (strippedEvent.length() < 2) {
            return false;
        }

        if (!Character.isLetter(strippedEvent.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(strippedEvent.charAt(1))) {
            return false;
        }

        // Some characters are definitely not single codon ranges
        if (strippedEvent.contains("*") || strippedEvent.contains("/") || strippedEvent.contains("fs")
                || strippedEvent.contains(",")) {
            return false;
        }

        // Some events contain multiple digit sequences.
        if (countDigitSequences(strippedEvent) > 1) {
            return false;
        }

        // Codon range should end with X or with a digit!
        return Character.isDigit(strippedEvent.charAt(strippedEvent.length() - 1)) || strippedEvent.endsWith("X");
    }

    @VisibleForTesting
    static int countDigitSequences(@NotNull String string) {
        int digitSequences = 0;

        boolean inDigitSequence = false;
        for (int i = 0; i < string.length(); i++) {
            if (Character.isDigit(string.charAt(i)) && !inDigitSequence) {
                inDigitSequence = true;
                digitSequences++;
            } else if (!Character.isDigit(string.charAt(i)) && inDigitSequence) {
                inDigitSequence = false;
            }
        }
        return digitSequences;
    }
}
