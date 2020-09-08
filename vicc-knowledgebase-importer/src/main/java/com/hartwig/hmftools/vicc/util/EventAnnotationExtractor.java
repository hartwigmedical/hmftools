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
            @Nullable String provenanceRule, @NotNull String proteinAnnotation) {
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

        if (DetermineHotspot.isResolvableProteinAnnotation(proteinAnnotation)) {
            return EventAnnotation.HOTSPOT;
        } else if (EventAnnotationExtractor.SIGNATURES.contains(feature)) {
            return EventAnnotation.SIGNATURE;
        } else if (DetermineCopyNumber.isAmplification(feature, biomarkerType)) {
            return EventAnnotation.AMPLIFICATION;
        } else if (DetermineCopyNumber.isDeletion(feature, biomarkerType)) {
            return EventAnnotation.DELETION;
        } else if (DetermineFusion.isFusion(feature, biomarkerType, provenanceRule, proteinAnnotation)) {
            return EventAnnotation.FUSION_PAIR;
        } else if (DetermineFusion.isFusionPromiscuous(feature, biomarkerType, provenanceRule, proteinAnnotation)) {
            return EventAnnotation.FUSION_PROMISCUOUS;
        } else if (!DetermineHotspot.isResolvableProteinAnnotation(proteinAnnotation)) {
            if (EventAnnotationExtractor.GENE_LEVEL.contains(biomarkerType) || EventAnnotationExtractor.GENE_LEVEL.contains(feature)
                    || EventAnnotationExtractor.GENE_LEVEL.contains(provenanceRule) || EventAnnotationExtractor.GENE_LEVEL.contains(
                    proteinAnnotation)) {
                return EventAnnotation.GENE_LEVEL;
            }
        } else if (EventAnnotationExtractor.GENE_EXON.contains(event) && !feature.toLowerCase().contains("deletion")) {
            return EventAnnotation.GENE_RANGE_EXON;
        } else if (EventAnnotationExtractor.GENE_MULTIPLE_CODONS.contains(biomarkerType) && proteinAnnotation.substring(
                proteinAnnotation.length() - 1).equals("X") && EventAnnotationExtractor.GENE_MULTIPLE_CODONS.contains(proteinAnnotation)) {
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
            LOGGER.warn("No event annotation extracted from event!");
            return EventAnnotation.UNKNOWN;
        }
        return EventAnnotation.UNKNOWN;

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
