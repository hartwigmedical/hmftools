package com.hartwig.hmftools.vicc.util;

import java.util.Set;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class EventAnnotationExtractor {

    private EventAnnotationExtractor() {
    }

    private static final Logger LOGGER = LogManager.getLogger(EventAnnotationExtractor.class);

    private static final String HGVS_RANGE_INDICATOR = "_";
    private static final String HGVS_DELETION = "del";
    private static final String HGVS_INSERTION = "ins";
    private static final String HGVS_DUPLICATION = "dup";

    private static final String HGVS_FRAMESHIFT_SUFFIX = "fs";
    private static final String HGVS_FRAMESHIFT_SUFFIX_WITH_STOP_GAINED = "fs*";

    private static final int MAX_INFRAME_BASE_LENGTH = 50;

    public static final Set<String> AMPLIFICATIONS = Sets.newHashSet("Amplification",
            "amplification",
            "AMPLIFICATION",
            "amp",
            "overexpression",
            "over exp",
            "amp over exp",
            "OVEREXPRESSION",
            "Overexpression");

    public static final Set<String> DELETIONS = Sets.newHashSet("Deletion",
            "deletion",
            "DELETION",
            "del",
            "undexpression",
            "dec exp",
            "UNDEREXPRESSION",
            "loss",
            "LOSS",
            "Copy Number Loss");

    public static final Set<String> SEARCH_FUSION_PAIRS = Sets.newHashSet("Fusion",
            "Disruptive Inframe Deletion",
            "Gene Fusion",
            "EGFR-KDD",
            "Transcript Regulatory Region Fusion",
            "FGFR3 - BAIAP2L1 Fusion",
            "p61BRAF-V600E");
    public static final Set<String> SEARCH_FUSION_PROMISCUOUS =
            Sets.newHashSet("REARRANGEMENT", "Fusions", "fusion", "rearrange", "Transcript Fusion", "FUSION", "nonsense", "FUSIONS");

    public static final Set<String> IGNORE = Sets.newHashSet("3' EXON DELETION");

    public static final Set<String> INTERNAL_FUSION =
            Sets.newHashSet("(Partial", "Exon Loss Variant", "Inframe Deletion", "is_deletion", "EGFRvIII", "EGFRvV", "EGFRvII", "ITD");

    public static final Set<String> GENE_LEVEL = Sets.newHashSet("Gain-of-function Mutations",
            "Gain-of-Function",
            "Oncogenic Mutations",
            "MUTATION",
            "act mut",
            "pos",
            "positive",
            "inact mut",
            "biallelic inactivation",
            "negative",
            "Loss Of Function Variant",
            "Loss Of Heterozygosity",
            "TRUNCATING MUTATION",
            "Truncating Mutations",
            "mutant",
            "mut",
            "gene_only",
            "ACTIVATING MUTATION",
            "DELETERIOUS MUTATION",
            "feature_truncation",
            "FRAMESHIFT TRUNCATION",
            "FRAMESHIFT MUTATION",
            "SPLICE VARIANT 7",
            "Splice",
            "DNMT3B7",
            "LCS6-variant",
            "AR-V7",
            "ARv567es");

    public static final Set<String> GENE_EXON = Sets.newHashSet("exon", "Exon Variant");
    public static final Set<String> GENE_MULTIPLE_CODONS =
            Sets.newHashSet("nonsense", "(V600)", "splice_region_variant", "Splice Donor Variant", "Inframe Deletion");

    public static final Set<String> SIGNATURES = Sets.newHashSet("Microsatellite Instability-High");

    @NotNull
    public static EventAnnotation toEventAnnotation(@NotNull String featureName, @Nullable String biomarkerType,
            @NotNull String provenanceRule, @NotNull String proteinAnnotation) {
        String feature = featureName;
        if (feature.contains(" ") && !feature.equals("Copy Number Loss")) {
            feature = feature.split(" ", 2)[1];
        }

        String event = Strings.EMPTY;
        if (feature.toLowerCase().contains("exon")) {
            event = "exon";
        } else if (biomarkerType != null) {
            if (biomarkerType.equals("Exon Variant")) {
                event = "exon";
            }
        }

        if (EventAnnotationExtractor.SIGNATURES.contains(feature)) {
            return EventAnnotation.SIGNATURE;
        } else if (EventAnnotationExtractor.AMPLIFICATIONS.contains(feature) || EventAnnotationExtractor.AMPLIFICATIONS.contains(
                biomarkerType)) {
            return EventAnnotation.AMPLIFICATION;
        } else if (EventAnnotationExtractor.DELETIONS.contains(feature) || EventAnnotationExtractor.DELETIONS.contains(biomarkerType)) {
            return EventAnnotation.DELETION;
        } else if (isFusion(feature, biomarkerType, provenanceRule, proteinAnnotation)) {
            return EventAnnotation.FUSION_PAIR;
        } else if (isFusionPromiscuous(feature, biomarkerType, provenanceRule, proteinAnnotation)) {
            return EventAnnotation.FUSION_PROMISCUOUS;
        } else if (isResolvableProteinAnnotation(proteinAnnotation)) {
            return EventAnnotation.HOTSPOT;
        } else if (!isResolvableProteinAnnotation(proteinAnnotation)) {
            if (EventAnnotationExtractor.GENE_LEVEL.contains(biomarkerType) || EventAnnotationExtractor.GENE_LEVEL.contains(feature)
                    || EventAnnotationExtractor.GENE_LEVEL.contains(provenanceRule) || EventAnnotationExtractor.GENE_LEVEL.contains(
                    proteinAnnotation)) {
                return EventAnnotation.GENE_LEVEL;
            } else {
                //            LOGGER.warn("No event annotation extracted from event!");
                return EventAnnotation.UNKNOWN;
            }
        } else if (EventAnnotationExtractor.GENE_EXON.contains(event) && !feature.toLowerCase().contains("deletion")) {
            return EventAnnotation.GENE_RANGE_EXON;
        } else if (EventAnnotationExtractor.GENE_MULTIPLE_CODONS.contains(biomarkerType) && proteinAnnotation
                .substring(proteinAnnotation.length() - 1)
                .equals("X") && EventAnnotationExtractor.GENE_MULTIPLE_CODONS.contains(proteinAnnotation)) {
            return EventAnnotation.GENE_RANGE_CODON;
        } else if (proteinAnnotation.length() >= 1 && isValidSingleCodonRange(proteinAnnotation)) {
            return EventAnnotation.GENE_RANGE_CODON;
        } else if (EventAnnotationExtractor.GENE_MULTIPLE_CODONS.contains(biomarkerType)) {
            return EventAnnotation.GENE_RANGE_CODON;
        } else if (feature.contains("DEL") && EventAnnotationExtractor.GENE_MULTIPLE_CODONS.contains(biomarkerType)) {
            return EventAnnotation.GENE_RANGE_CODON;
        } else if (proteinAnnotation.contains("del") && proteinAnnotation.contains("_")) {
            return EventAnnotation.GENE_RANGE_CODON;
        } else {
            //            LOGGER.warn("No event annotation extracted from event!");
            return EventAnnotation.UNKNOWN;
        }

    }

    private static boolean isFusion(@NotNull String feature, @Nullable String biomarkerType, @NotNull String provenanceRule,
            @NotNull String proteinAnnotation) {
        FusionEvent eventKeyFusion = extractKeyFusion(feature, biomarkerType, provenanceRule, proteinAnnotation);
        if (eventKeyFusion.equals(FusionEvent.FUSION_PAIR)) {
            return true;
        } else {
            return false;
        }
    }

    private static boolean isFusionPromiscuous(@NotNull String feature, @Nullable String biomarkerType, @NotNull String provenanceRule,
            @NotNull String proteinAnnotation) {
        FusionEvent eventKeyFusion = extractKeyFusion(feature, biomarkerType, provenanceRule, proteinAnnotation);
        if (eventKeyFusion.equals(FusionEvent.FUSION_PROMISCUOUS)) {
            return true;
        } else {
            return false;
        }
    }

    @NotNull
    private static FusionEvent extractKeyFusion(@NotNull String feature, @Nullable String biomarkerType, @NotNull String provenanceRule,
            @NotNull String proteinAnnotation) {
        if (!EventAnnotationExtractor.IGNORE.contains(feature)) { // Extract internal fusion
            if (EventAnnotationExtractor.INTERNAL_FUSION.contains(proteinAnnotation)) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("DELETION") && EventAnnotationExtractor.INTERNAL_FUSION.contains(biomarkerType)) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("DELETION") && EventAnnotationExtractor.INTERNAL_FUSION.contains(provenanceRule)) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("Exon") && feature.contains("deletion")) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("Exon") && feature.contains("deletion/insertion")) {
                return FusionEvent.FUSION_PAIR;
            } else if (feature.contains("Exon") && feature.contains("insertions/deletions")) {
                return FusionEvent.FUSION_PAIR;
            } else if (EventAnnotationExtractor.SEARCH_FUSION_PAIRS.contains(proteinAnnotation)) {
                return FusionEvent.FUSION_PAIR;
            } else if (EventAnnotationExtractor.SEARCH_FUSION_PAIRS.contains(feature)) {
                return FusionEvent.FUSION_PAIR;
            } else if (EventAnnotationExtractor.SEARCH_FUSION_PROMISCUOUS.contains(proteinAnnotation)) {
                if (feature.contains("-")) {
                    return FusionEvent.FUSION_PAIR;
                } else {
                    return FusionEvent.FUSION_PROMISCUOUS;
                }
            } else if (biomarkerType != null) {
                if (EventAnnotationExtractor.SEARCH_FUSION_PROMISCUOUS.contains(biomarkerType)) {
                    if (feature.contains("-")) {
                        return FusionEvent.FUSION_PAIR;
                    } else {
                        return FusionEvent.FUSION_PROMISCUOUS;
                    }
                }
                if (EventAnnotationExtractor.SEARCH_FUSION_PAIRS.contains(biomarkerType)) {
                    return FusionEvent.FUSION_PAIR;
                }
            } else if (feature.toLowerCase().contains("exon") && feature.toLowerCase().contains("deletion") && feature.contains("-")
                    || feature.contains("&")) {// Extract internal fusion
                return FusionEvent.FUSION_PAIR;
            } else {
                return FusionEvent.UNKNOWN;
            }
        }

        //TODO: check why this is needed??
        return FusionEvent.UNKNOWN;
    }

    public static boolean isResolvableProteinAnnotation(@NotNull String proteinAnnotation) {
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
            LOGGER.warn("Could not determine whether protein annotation is resolvable due to '{}'", exception.getMessage(), exception);
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

            long start = Long.parseLong(annotationStartPart.substring(1));
            long end = Long.parseLong(annotationEndPart.substring(1, indexOfEvent));
            return 3 * (1 + end - start) <= MAX_INFRAME_BASE_LENGTH;
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

    private static boolean isInteger(@NotNull String value) {
        try {
            Integer.parseInt(value);
            return true;
        } catch (NumberFormatException exp) {
            return false;
        }
    }

    private static boolean isValidSingleCodonRange(@NotNull String feature) {

        // Features are expected to look something like V600 (1 char - N digits)
        if (feature.length() < 3) {
            return false;
        }

        if (!Character.isLetter(feature.charAt(0))) {
            return false;
        }

        if (!Character.isDigit(feature.charAt(1))) {
            return false;
        }

        if (feature.contains("*")) {
            return false;
        }

        if (feature.contains("/")) {
            return false;
        }

        if (feature.contains("fs")) {
            return false;
        }

        return Character.isDigit(feature.substring(feature.length() - 1).charAt(0));
    }

}
