package com.hartwig.hmftools.common.serve.classification.matchers;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;

class GeneRangeCodonMatcher implements EventMatcher {

    @NotNull
    private final EventPreprocessor proteinAnnotationExtractor;

    GeneRangeCodonMatcher(@NotNull final EventPreprocessor proteinAnnotationExtractor) {
        this.proteinAnnotationExtractor = proteinAnnotationExtractor;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        String processedEvent = proteinAnnotationExtractor.apply(event);

        // Feature codon ranges occasionally come with parentheses
        String strippedEvent = processedEvent.replace("(", "").replace(")", "");

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
        if (strippedEvent.contains("*") || strippedEvent.contains("/") || strippedEvent.contains("fs") || strippedEvent.contains(",")) {
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
