package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

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

    private static final String EXON_KEYWORD = "exon";
    private static final Set<String> EXON_RANGE_EXACT_TERMS = Sets.newHashSet("RARE EX 18-21 MUT");

    private static final Set<String> EXON_RANGE_KEYWORDS =
            Sets.newHashSet("deletion", "insertion", "proximal", "mutation", "splice site insertion", "frameshift");

    private GeneRangeClassifier() {
    }

    public static boolean isGeneLevelEvent(@NotNull String gene, @NotNull String event) {
        if (CombinedClassifier.isCombinedEvent(gene, event) || ComplexClassifier.isComplexEvent(gene, event)) {
            return false;
        }

        if (event.toLowerCase().contains(EXON_KEYWORD)) {
            return false;
        }

        for (String keyword : GENERIC_GENE_LEVEL_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : INACTIVATING_GENE_LEVEL_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        for (String keyword : ACTIVATING_GENE_LEVEL_KEYWORDS) {
            if (event.contains(keyword)) {
                return true;
            }
        }

        return event.trim().equals(gene);
    }

    public static boolean isGeneRangeExonEvent(@NotNull String gene, @NotNull String event) {
        if (CombinedClassifier.isFusionPairAndGeneRangeExon(gene, event) || CombinedClassifier.isCombinedEvent(gene, event)
                || ComplexClassifier.isComplexEvent(gene, event)) {
            return false;
        }

        if (EXON_RANGE_EXACT_TERMS.contains(event)) {
            return true;
        } else {
            String lowerCaseEvent = event.toLowerCase();
            if (lowerCaseEvent.contains(EXON_KEYWORD)) {
                for (String keyword : EXON_RANGE_KEYWORDS) {
                    if (lowerCaseEvent.contains(keyword)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    public static boolean isGeneRangeCodonEvent(@NotNull String event) {
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
