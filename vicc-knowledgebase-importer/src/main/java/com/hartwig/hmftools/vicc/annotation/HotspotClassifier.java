package com.hartwig.hmftools.vicc.annotation;

import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_DELETION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_DUPLICATION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_INSERTION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_RANGE_INDICATOR;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_START_LOST;

import java.util.List;

import com.hartwig.hmftools.common.serve.classification.CompositeEventMatcher;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HotspotClassifier implements EventMatcher {

    private static final int MAX_INFRAME_BASE_LENGTH = 50;

    @NotNull
    public static EventMatcher create(@NotNull List<EventMatcher> noMatchEventMatchers, @NotNull EventPreprocessor preprocessor) {
        return new CompositeEventMatcher(noMatchEventMatchers, new HotspotClassifier(preprocessor));
    }

    public static boolean isProteinAnnotation(@NotNull String event) {
        // Do note: No pre-processing is done on this event!
        return isProteinAnnotationWithMaxLength(event, null);
    }

    @NotNull
    private final EventPreprocessor preprocessor;

    private HotspotClassifier(@NotNull final EventPreprocessor preprocessor) {
        this.preprocessor = preprocessor;
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        String processedEvent = preprocessor.apply(event);

        if (isProteinAnnotationWithMaxLength(processedEvent, MAX_INFRAME_BASE_LENGTH)) {
            return !isHotspotOnFusionGene(event);
        }

        return false;
    }

    private static boolean isProteinAnnotationWithMaxLength(@NotNull String event, @Nullable Integer maxLength) {
        if (isValidFrameshift(event)) {
            return true;
        } else if (event.contains(HGVS_RANGE_INDICATOR)) {
            int mutationLength = extractRangeLength(event);
            return mutationLength > 0 && (maxLength == null || mutationLength <= maxLength);
        } else if (event.contains(HGVS_DELETION + HGVS_INSERTION)) {
            int mutationLength = extractComplexDeletionInsertionLength(event);
            return mutationLength > 0 && (maxLength == null || mutationLength <= maxLength);
        } else if (event.startsWith(HGVS_START_LOST)) {
            return true;
        } else {
            return isValidSingleCodonMutation(event);
        }
    }

    private static boolean isValidFrameshift(@NotNull String event) {
        if (event.endsWith(HGVS_FRAMESHIFT_SUFFIX) || event.endsWith(HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED)) {
            int frameshiftPosition = event.indexOf(HGVS_FRAMESHIFT_SUFFIX);
            if (frameshiftPosition > 1) {
                return isInteger(event.substring(frameshiftPosition - 1, frameshiftPosition));
            }
        }

        return false;
    }

    private static int extractRangeLength(@NotNull String event) {
        assert event.contains(HGVS_RANGE_INDICATOR);

        // Features could be ranges such as E102_I103del. We whitelist specific feature types when analyzing a range.
        String[] parts = event.split(HGVS_RANGE_INDICATOR);
        String startPart = parts[0];
        String endPart = parts[1];
        if (endPart.contains(HGVS_INSERTION) || endPart.contains(HGVS_DUPLICATION) || endPart.contains(HGVS_DELETION)) {
            int indexOfEvent;
            // Keep in mind that 'del' always comes prior to 'ins' in situations of complex inframes.
            if (endPart.contains(HGVS_DELETION)) {
                indexOfEvent = endPart.indexOf(HGVS_DELETION);
            } else if (endPart.contains(HGVS_DUPLICATION)) {
                indexOfEvent = endPart.indexOf(HGVS_DUPLICATION);
            } else {
                indexOfEvent = endPart.indexOf(HGVS_INSERTION);
            }

            String startRange = startPart.substring(1);
            String endRange = endPart.substring(1, indexOfEvent);

            if (isLong(startRange) && isLong(endRange)) {
                return (int) (3 * (1 + Long.parseLong(endRange) - Long.parseLong(startRange)));
            } else {
                return -1;
            }
        } else {
            return -1;
        }
    }

    private static int extractComplexDeletionInsertionLength(@NotNull String event) {
        assert event.contains(HGVS_DELETION + HGVS_INSERTION);

        String[] parts = event.split(HGVS_DELETION + HGVS_INSERTION);

        // Format is expected to be something like D770delinsGY
        if (isInteger(parts[0].substring(1))) {
            return 3 * parts[1].length();
        } else {
            return -1;
        }
    }

    private static boolean isValidSingleCodonMutation(@NotNull String event) {
        // Single codon mutations are expected to look something like V600E (1 char - N digits - M chars (1 char, or "del" or "dup"))
        if (event.length() < 3) {
            return false;
        }

        if (!Character.isLetter(event.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(event.charAt(1))) {
            return false;
        }

        boolean haveObservedNonDigit = !Character.isDigit(event.charAt(2));
        int firstNotDigit = haveObservedNonDigit ? 2 : -1;
        for (int i = 3; i < event.length(); i++) {
            char charToEvaluate = event.charAt(i);
            if (haveObservedNonDigit && Character.isDigit(charToEvaluate)) {
                return false;
            }
            boolean isDigit = Character.isDigit(charToEvaluate);
            if (!isDigit && firstNotDigit == -1) {
                firstNotDigit = i;
            }

            haveObservedNonDigit = haveObservedNonDigit || !isDigit;
        }

        if (!haveObservedNonDigit) {
            return false;
        }

        String newAminoAcid = event.substring(firstNotDigit);
        // X is a wildcard which should become a codon range rather than a single codon mutation.
        return !newAminoAcid.equals("X") && (newAminoAcid.length() == 1 || newAminoAcid.equals(HGVS_DELETION) || newAminoAcid.equals(
                HGVS_DUPLICATION));
    }

    private static boolean isHotspotOnFusionGene(@NotNull String event) {
        String trimmedEvent = event.trim();
        if (trimmedEvent.contains(" ")) {
            String[] parts = trimmedEvent.split(" ");
            return FusionPairClassifier.isFusionPair(parts[0]);
        }
        return false;
    }

    private static boolean isLong(@NotNull String value) {
        try {
            Long.parseLong(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
    }

    private static boolean isInteger(@NotNull String value) {
        try {
            Integer.parseInt(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
    }
}
