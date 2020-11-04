package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class CopyNumberClassifier {

    private static final Set<String> AMPLIFICATION_KEYWORDS = Sets.newHashSet("Amplification",
            "amplification",
            "AMPLIFICATION",
            "amp",
            "overexpression",
            "over exp",
            "amp over exp",
            "OVEREXPRESSION",
            "Overexpression");

    private static final Set<String> DELETION_KEYWORDS = Sets.newHashSet("Deletion",
            "deletion",
            "DELETION",
            "del",
            "undexpression",
            "dec exp",
            "UNDEREXPRESSION",
            "loss",
            "LOSS",
            "Copy Number Loss");

    private CopyNumberClassifier() {
    }

    public static boolean isAmplification(@NotNull String featureName, @Nullable String biomarkerType) {
        return extractAmplificationDeletionType(featureName, biomarkerType) == CopyNumberType.AMPLIFICATION;
    }

    public static boolean isDeletion(@NotNull String featureName, @Nullable String biomarkerType) {
        if (!featureName.toLowerCase().contains("exon")) {
            return extractAmplificationDeletionType(featureName, biomarkerType) == CopyNumberType.DELETION;
        }

        return false;
    }

    @Nullable
    private static CopyNumberType extractAmplificationDeletionType(@NotNull String featureName, @Nullable String biomarkerType) {
        if (featureName.contains(" ")) {
            // Never consider something an amp or del when the first part has a dash.
            if (featureName.split(" ")[0].contains("-")) {
                return null;
            }

        }
        for (String keyword : AMPLIFICATION_KEYWORDS) {
            if (featureName.contains(keyword)) {
                return CopyNumberType.AMPLIFICATION;
            }
        }

        for (String keyword : DELETION_KEYWORDS) {
            if (featureName.contains(keyword)) {
                return CopyNumberType.DELETION;
            }
        }

        return null;
    }

    private enum CopyNumberType {
        AMPLIFICATION,
        DELETION
    }
}
