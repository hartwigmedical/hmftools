package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FeatureTypeExtractor {

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
            "fusion",
            "EGFR-KDD",
            "Transcript Regulatory Region Fusion",
            "FGFR3 - BAIAP2L1 Fusion",
            "p61BRAF-V600E");
    public static final Set<String> SEARCH_FUSION_PROMISCUOUS =
            Sets.newHashSet("REARRANGEMENT", "Fusions", "fusion", "rearrange", "Transcript Fusion", "FUSION", "FUSIONS");

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
            Sets.newHashSet("(V600)", "splice_region_variant", "Splice Donor Variant", "Inframe Deletion");

    public static final Set<String> SIGNATURES = Sets.newHashSet("Microsatellite Instability-High");

    private FeatureTypeExtractor() {
    }

    @NotNull
    public static FeatureType extractType(@NotNull Feature feature) {
        return extractType(feature.name(),
                feature.biomarkerType(),
                feature.provenanceRule(),
                ProteinAnnotationExtractor.extractProteinAnnotation(feature));
    }

    @NotNull
    public static FeatureType extractType(@NotNull String featureName, @Nullable String biomarkerType, @Nullable String provenanceRule,
            @NotNull String proteinAnnotation) {
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

        if (biomarkerType != null && provenanceRule != null) {
            if (featureName.contains("+") && (biomarkerType.equals("amp") && provenanceRule.contains("is_fusion_acceptor") || provenanceRule
                    .contains("is_fusion_donor"))) {
                return FeatureType.COMBINED;
            } else if (featureName.contains("insertion")) {
                int countInsertion = featureName.split("insertion").length -1;
                if (countInsertion > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("deletion")) {
                int countDeletion = featureName.split("deletion").length -1;
                if (countDeletion > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("frameshift")) {
                int countFrameshift = featureName.split("frameshift").length -1;
                if (countFrameshift > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("insertions") && featureName.contains("deletion")) {
                int countCombined = (featureName.split("insertion").length -1) + (featureName.split("deletion").length -1);
                if (countCombined > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("splice")) {
                int countSplice = featureName.split("splice").length-1;
                if (countSplice >1) {
                    return FeatureType.COMBINED;
                }
            }
        }
        if (DetermineHotspot.isHotspot(proteinAnnotation)) {
            return FeatureType.HOTSPOT;
        } else if (FeatureTypeExtractor.SIGNATURES.contains(feature)) {
            return FeatureType.SIGNATURE;
        } else if (DetermineCopyNumber.isAmplification(feature, biomarkerType)) {
            return FeatureType.AMPLIFICATION;
        } else if (DetermineCopyNumber.isDeletion(feature, biomarkerType) && !feature.toLowerCase().contains("exon")) {
            return FeatureType.DELETION;
        } else if (DetermineFusion.isFusion(feature, biomarkerType, provenanceRule, proteinAnnotation)) {
            return FeatureType.FUSION_PAIR;
        } else if (DetermineFusion.isFusionPromiscuous(feature, biomarkerType, provenanceRule, proteinAnnotation)) {
            return FeatureType.FUSION_PROMISCUOUS;
        }  else if (FeatureTypeExtractor.GENE_EXON.contains(event) && !feature.toLowerCase().contains("deletion")) {
            return FeatureType.GENE_RANGE_EXON;
        } else if (FeatureTypeExtractor.GENE_MULTIPLE_CODONS.contains(biomarkerType) && proteinAnnotation.substring(
                proteinAnnotation.length() - 1).equals("X") && FeatureTypeExtractor.GENE_MULTIPLE_CODONS.contains(proteinAnnotation)) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (proteinAnnotation.length() >= 1 && isValidSingleCodonRange(proteinAnnotation)) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (FeatureTypeExtractor.GENE_MULTIPLE_CODONS.contains(biomarkerType)) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (feature.contains("DEL") && FeatureTypeExtractor.GENE_MULTIPLE_CODONS.contains(biomarkerType)) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (proteinAnnotation.contains("del") && proteinAnnotation.contains("_")) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (!DetermineHotspot.isHotspot(proteinAnnotation)) {
            if (FeatureTypeExtractor.GENE_LEVEL.contains(biomarkerType) || FeatureTypeExtractor.GENE_LEVEL.contains(feature)
                    || FeatureTypeExtractor.GENE_LEVEL.contains(provenanceRule) || FeatureTypeExtractor.GENE_LEVEL.contains(
                    proteinAnnotation)) {
                return FeatureType.GENE_LEVEL;
            }
        }

        return FeatureType.UNKNOWN;
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
