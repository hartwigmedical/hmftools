package com.hartwig.hmftools.vicc.annotation;

import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_DELETION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_DUPLICATION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_INSERTION;
import static com.hartwig.hmftools.common.variant.hgvs.HgvsConstants.HGVS_RANGE_INDICATOR;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class DetermineHotspot {

    private DetermineHotspot() {
    }

    private static final Logger LOGGER = LogManager.getLogger(DetermineHotspot.class);

    private static final int MAX_INFRAME_BASE_LENGTH = 50;

    public static boolean isHotspot(@NotNull String proteinAnnotation) {
        try {
            if (isFrameshift(proteinAnnotation)) {
                return isValidFrameshift(proteinAnnotation);
            } else if (proteinAnnotation.contains(HGVS_RANGE_INDICATOR)) {
                return isValidRangeMutation(proteinAnnotation);
            } else if (proteinAnnotation.contains(HGVS_DELETION + HGVS_INSERTION)) {
                return isValidComplexDeletionInsertion(proteinAnnotation);
            } else if (proteinAnnotation.startsWith("*")) {
                return true;
            } else {
                return isValidSingleCodonMutation(proteinAnnotation);
            }
        } catch (Exception exception) {
            LOGGER.warn("Could not determine whether protein annotation '{}' is resolvable due to '{}'",
                    proteinAnnotation,
                    exception.getMessage(),
                    exception);
            return false;
        }
    }

    private static boolean isFrameshift(@NotNull String proteinAnnotation) {
        return proteinAnnotation.endsWith(HGVS_FRAMESHIFT_SUFFIX) || proteinAnnotation.endsWith(HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED);
    }

    private static boolean isValidFrameshift(@NotNull String proteinAnnotation) {
        int frameshiftPosition = proteinAnnotation.indexOf(HGVS_FRAMESHIFT_SUFFIX);
        if (frameshiftPosition > 1) {
            return isInteger(proteinAnnotation.substring(frameshiftPosition - 1, frameshiftPosition));
        }

        return false;
    }

    private static boolean isValidRangeMutation(@NotNull String proteinAnnotation) {
        assert proteinAnnotation.contains(HGVS_RANGE_INDICATOR);

        // Features could be ranges such as E102_I103del. We whitelist specific feature types when analyzing a range.
        String[] annotationParts = proteinAnnotation.split(HGVS_RANGE_INDICATOR);
        String annotationStartPart = annotationParts[0];
        String annotationEndPart = annotationParts[1];
        if (annotationEndPart.contains(HGVS_INSERTION) || annotationEndPart.contains(HGVS_DUPLICATION) || annotationEndPart.contains(
                HGVS_DELETION)) {
            int indexOfEvent;
            // Keep in mind that 'del' always comes prior to 'ins' in situations of complex inframes.
            if (annotationEndPart.contains(HGVS_DELETION)) {
                indexOfEvent = annotationEndPart.indexOf(HGVS_DELETION);
            } else if (annotationEndPart.contains(HGVS_DUPLICATION)) {
                indexOfEvent = annotationEndPart.indexOf(HGVS_DUPLICATION);
            } else {
                indexOfEvent = annotationEndPart.indexOf(HGVS_INSERTION);
            }

            String startRange = annotationStartPart.substring(1);
            String endRange = annotationEndPart.substring(1, indexOfEvent);

            if (isLong(startRange) && isLong(endRange)) {
                return 3 * (1 + Long.parseLong(endRange) - Long.parseLong(startRange)) <= MAX_INFRAME_BASE_LENGTH;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    private static boolean isValidComplexDeletionInsertion(@NotNull String proteinAnnotation) {
        String[] annotationParts = proteinAnnotation.split(HGVS_DELETION + HGVS_INSERTION);

        return isInteger(annotationParts[0].substring(1)) && (3 * annotationParts[1].length()) <= MAX_INFRAME_BASE_LENGTH;
    }

    private static boolean isValidSingleCodonMutation(@NotNull String proteinAnnotation) {
        if (proteinAnnotation.contains(HGVS_INSERTION)) {
            // Insertions are only allowed in a range, since we need to know where to insert the sequence exactly.
            return false;
        }

        // Features are expected to look something like V600E (1 char - N digits - M chars)
        if (proteinAnnotation.length() < 3) {
            return false;
        }

        if (!Character.isLetter(proteinAnnotation.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(proteinAnnotation.charAt(1))) {
            return false;
        }

        boolean haveObservedNonDigit = !Character.isDigit(proteinAnnotation.charAt(2));
        int firstNotDigit = haveObservedNonDigit ? 2 : -1;
        for (int i = 3; i < proteinAnnotation.length(); i++) {
            char charToEvaluate = proteinAnnotation.charAt(i);
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

        String newAminoAcid = proteinAnnotation.substring(firstNotDigit);
        // X is a wildcard that we don't support, and "/" indicates logical OR that we don't support.
        return !newAminoAcid.equals("X") && !newAminoAcid.contains("/");
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
