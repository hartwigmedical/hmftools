package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class CopyNumberClassifier {

    private static final Set<String> AMPLIFICATIONS = Sets.newHashSet("Amplification",
            "amplification",
            "AMPLIFICATION",
            "amp",
            "overexpression",
            "over exp",
            "amp over exp",
            "OVEREXPRESSION",
            "Overexpression");

    private static final Set<String> DELETIONS = Sets.newHashSet("Deletion",
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
        String correctedFeatureName;
        if (featureName.contains(" ") && !featureName.equals("Copy Number Loss")) {
            correctedFeatureName = featureName.split(" ", 2)[1];
        } else {
            correctedFeatureName = featureName;
        }

        if (AMPLIFICATIONS.contains(correctedFeatureName) || AMPLIFICATIONS.contains(biomarkerType)) {
            return CopyNumberType.AMPLIFICATION;
        } else if (DELETIONS.contains(correctedFeatureName) || DELETIONS.contains(biomarkerType)) {
            return CopyNumberType.DELETION;
        }

        return null;
    }

    private enum CopyNumberType {
        AMPLIFICATION,
        DELETION
    }
}
