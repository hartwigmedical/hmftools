package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GeneRangeClassifier {

    public static final Set<String> DETAILED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO = Sets.newHashSet("MUTATION",
            "mutant",
            "mut",
            "TRUNCATING MUTATION",
            "Truncating Mutations",
            "feature_truncation",
            "FRAMESHIFT TRUNCATION",
            "FRAMESHIFT MUTATION",
            "ALTERATION");

    public static final Set<String> DETAILED_GENE_LEVEL_INFO_WITH_TSG = Sets.newHashSet("inact mut",
            "biallelic inactivation",
            "Loss Of Function Variant",
            "Loss Of Heterozygosity",
            "DELETERIOUS MUTATION",
            "negative",
            "BIALLELIC INACTIVATION",
            "LOSS-OF-FUNCTION");

    public static final Set<String> DETAILED_GENE_LEVEL_INFO_WITH_ONCO = Sets.newHashSet("Gain-of-function Mutations",
            "Gain-of-Function",
            "act mut",
            "ACTIVATING MUTATION", "Oncogenic Mutations", "pos", "positive", "oncogenic mutation");

    public static final String GENE_LEVEL = "gene_only";

    private GeneRangeClassifier() {
    }

    public static boolean isGeneLevelEvent(@NotNull String featureName, @Nullable String provenanceRule) {
        String[] parts = featureName.split(" ");
        if (parts.length < 2) {
            return false;
        }
        String eventDescription = parts[1].trim();

        return DETAILED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription) || DETAILED_GENE_LEVEL_INFO_WITH_TSG.contains(
                eventDescription) || DETAILED_GENE_LEVEL_INFO_WITH_ONCO.contains(eventDescription) || GENE_LEVEL.equals(provenanceRule);
    }

    public static boolean isGeneRangeExonEvent(@NotNull String featureName, @NotNull String gene) {
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
