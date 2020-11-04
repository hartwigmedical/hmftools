package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GeneRangeClassifier {

    public static final Set<String> GENERIC_GENE_LEVEL_KEYWORDS = Sets.newHashSet("MUTATION",
            "mutant",
            "mut",
            "TRUNCATING MUTATION",
            "Truncating Mutations",
            "feature_truncation",
            "FRAMESHIFT TRUNCATION",
            "FRAMESHIFT MUTATION",
            "ALTERATION");

    public static final Set<String> INACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet("inact mut",
            "biallelic inactivation",
            "Loss Of Function Variant",
            "Loss Of Heterozygosity",
            "DELETERIOUS MUTATION",
            "negative",
            "BIALLELIC INACTIVATION",
            "LOSS-OF-FUNCTION");

    public static final Set<String> ACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet("Gain-of-function Mutations",
            "Gain-of-Function",
            "act mut",
            "ACTIVATING MUTATION",
            "Oncogenic Mutations",
            "pos",
            "positive",
            "oncogenic mutation");

    public static final String GENE_ONLY = "gene_only";

    private GeneRangeClassifier() {
    }

    public static boolean isGeneLevelEvent(@NotNull String featureName, @Nullable String provenanceRule) {
        for (String keyword : GENERIC_GENE_LEVEL_KEYWORDS) {
            if (featureName.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : INACTIVATING_GENE_LEVEL_KEYWORDS) {
            if (featureName.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : ACTIVATING_GENE_LEVEL_KEYWORDS) {
            if (featureName.contains(keyword)) {
                return true;
            }
        }

        return provenanceRule != null && provenanceRule.equals(GENE_ONLY);
    }

    public static boolean isGeneRangeExonEvent(@NotNull String featureName, @Nullable String gene) {
        if (CombinedClassifier.isFusionPairAndGeneRangeExon(featureName, gene) || CombinedClassifier.isCombinedEvent(featureName)) {
            return false;
        }

        String lowerCaseFeatureName = featureName.toLowerCase();
        if (lowerCaseFeatureName.contains("exon")) {
            return lowerCaseFeatureName.contains("deletion") || lowerCaseFeatureName.contains("insertion") || lowerCaseFeatureName.contains(
                    "proximal") || lowerCaseFeatureName.contains("mutation") || lowerCaseFeatureName.contains("splice site insertion")
                    || lowerCaseFeatureName.contains("frameshift");
        }
        return false;
    }

    public static boolean isGeneRangeCodonEvent(@NotNull String featureName) {
        if (featureName.length() > 1 && featureName.substring(featureName.length() - 1).equals("X")) {
            return true;
        } else {
            return featureName.length() >= 1 && isValidSingleCodonRange(featureName);
        }
    }

    private static boolean isValidSingleCodonRange(@NotNull String featureName) {
        // Features are expected to look something like V600 (1 char - N digits)
        if (featureName.length() < 3) {
            return false;
        }

        if (!Character.isLetter(featureName.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(featureName.charAt(1))) {
            return false;
        }

        if (featureName.contains("*")) {
            return false;
        }

        if (featureName.contains("/")) {
            return false;
        }

        if (featureName.contains("fs")) {
            return false;
        }

        return Character.isDigit(featureName.substring(featureName.length() - 1).charAt(0));
    }
}
