package com.hartwig.hmftools.vicc.util;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DetermineCopyNumber {

    private DetermineCopyNumber() {
    }

    public static boolean isAmplification(@NotNull String feature, @Nullable String biomarkertype) {
        String eventKeyAmplification = extractKeyAmplificationDeletion(feature, biomarkertype);
        if (eventKeyAmplification.equals("Amplification")) {
            return true;
        } else {
            return false;
        }
    }

    public static boolean isDeletion(@NotNull String feature, @Nullable String biomarkertype) {
        String eventKeyDeletion = extractKeyAmplificationDeletion(feature, biomarkertype);

        if (eventKeyDeletion.equals("Deletion")) {
            return true;
        } else {
            return false;
        }
    }

    private static String extractKeyAmplificationDeletion(@NotNull String feature, @Nullable String biomarkertype) {
        //TODO: fix combi events
        String featureName = feature;
        if (featureName.contains(" ") && !featureName.equals("Copy Number Loss")) {
            featureName = featureName.split(" ", 2)[1];
        }

        if (EventAnnotationExtractor.AMPLIFICATIONS.contains(featureName)
                || EventAnnotationExtractor.AMPLIFICATIONS.contains(biomarkertype)) {
            return "Amplification";
        } else if (EventAnnotationExtractor.DELETIONS.contains(featureName) || EventAnnotationExtractor.DELETIONS.contains(biomarkertype)) {
            return "Deletion";
        } else {
            return Strings.EMPTY;
        }
    }
}
