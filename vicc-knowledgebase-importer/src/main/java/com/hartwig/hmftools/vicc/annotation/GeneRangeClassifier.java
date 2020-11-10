package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
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

    private GeneRangeClassifier() {
    }

    public static boolean isGeneLevelEvent(@NotNull String featureName, @Nullable String gene) {
        if (CombinedClassifier.isCombinedEvent(featureName, gene) || ComplexClassifier.isComplexEvent(featureName, gene)) {
            return false;
        }

        if (featureName.toLowerCase().contains("exon")) {
            return false;
        }

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

        return featureName.trim().equals(gene);
    }

    public static boolean isGeneRangeExonEvent(@NotNull String featureName, @Nullable String gene) {
        if (CombinedClassifier.isFusionPairAndGeneRangeExon(featureName, gene) || CombinedClassifier.isCombinedEvent(featureName, gene)) {
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
        String proteinAnnotation = HotspotClassifier.extractProteinAnnotation(featureName);

        return isValidSingleCodonRange(proteinAnnotation);
    }

    private static boolean isValidSingleCodonRange(@NotNull String featureName) {
        // Feature codon ranges occasionally come with parentheses
        String strippedFeature = featureName.replace("(", "").replace(")", "");

        // Features are expected to look something like V600 (1 char - N digits)
        if (strippedFeature.length() < 2) {
            return false;
        }

        if (!Character.isLetter(strippedFeature.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(strippedFeature.charAt(1))) {
            return false;
        }

        // Some characters are definitely not single codon ranges
        if (strippedFeature.contains("*") || strippedFeature.contains("/") || strippedFeature.contains("fs")
                || strippedFeature.contains(",")) {
            return false;
        }

        // Some features contain multiple digit sequences.
        if (countDigitSequences(strippedFeature) > 1) {
            return false;
        }

        // Codon range should end with X or with a digit!
        return Character.isDigit(strippedFeature.charAt(strippedFeature.length() - 1)) || strippedFeature.endsWith("X");
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
