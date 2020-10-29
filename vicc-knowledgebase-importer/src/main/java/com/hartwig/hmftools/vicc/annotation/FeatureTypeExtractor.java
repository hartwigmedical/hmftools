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
            "FLT3-ITD");
    public static final Set<String> SEARCH_FUSION_PROMISCUOUS =
            Sets.newHashSet("REARRANGEMENT", "Fusions", "fusion", "rearrange", "Transcript Fusion", "FUSION", "FUSIONS");

    public static final Set<String> IGNORE = Sets.newHashSet("3' EXON DELETION");

    public static final Set<String> INTERNAL_FUSION = Sets.newHashSet("is_deletion", "EGFRvIII", "EGFRvV", "EGFRvII", "ITD");

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
            "ACTIVATING MUTATION",
            "Oncogenic Mutations",
            "pos",
            "positive",
            "oncogenic mutation");

    public static final Set<String> GENE_LEVEL = Sets.newHashSet("gene_only");

    public static final Set<String> GENE_EXON = Sets.newHashSet("exon");

    public static final Set<String> SIGNATURES = Sets.newHashSet("Microsatellite Instability-High");

    //TODO: check if more EXON_DEL_DUP fusions need to be added
    public static final Set<String> FUSION_PAIR_AND_EXON_RANGE =
            Sets.newHashSet("KIT EXON 11 MUTATION", "KIT Exon 11 mutations", "KIT Exon 11 deletions", "MET EXON 14 SKIPPING MUTATION");

    private FeatureTypeExtractor() {
    }

    @NotNull
    public static FeatureType extractType(@NotNull Feature feature) {
        return extractType(feature.name(),
                feature.biomarkerType(),
                feature.provenanceRule(),
                feature.description());
    }

    @NotNull
    public static FeatureType extractType(@NotNull String featureName, @Nullable String biomarkerType, @Nullable String provenanceRule,
            @Nullable String featureDescription) {
        String proteinAnnotation = DetermineHotspot.extractProteinAnnotation(featureName);
        String event = Strings.EMPTY;
        if (featureName.toLowerCase().contains("exon")) {
            event = "exon";
        }

        if (biomarkerType != null && provenanceRule != null) {
            if (featureName.contains("+") && !featureName.toLowerCase().contains("c.") && !featureName.contains(">")) {
                return FeatureType.COMBINED;
            } else if (featureName.contains("insertion")) {
                int countInsertion = featureName.split("insertion").length - 1;
                if (countInsertion > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("deletion")) {
                int countDeletion = featureName.split("deletion").length - 1;
                if (countDeletion > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("frameshift")) {
                int countFrameshift = featureName.split("frameshift").length - 1;
                if (countFrameshift > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("insertions") && featureName.contains("deletion")) {
                int countCombined = (featureName.split("insertion").length - 1) + (featureName.split("deletion").length - 1);
                if (countCombined > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (featureName.contains("splice")) {
                int countSplice = featureName.split("splice").length - 1;
                if (countSplice > 1) {
                    return FeatureType.COMBINED;
                }
            } else if (proteinAnnotation.equals("p61BRAF-V600E")) {
                return FeatureType.COMBINED;
            }
        }
        if (DetermineHotspot.isHotspot(featureName)) {
            return FeatureType.HOTSPOT;
        } else if (FeatureTypeExtractor.FUSION_PAIR_AND_EXON_RANGE.contains(featureDescription)) {
            return FeatureType.FUSION_PAIR_AND_GENE_RANGE_EXON;
        } else if (FeatureTypeExtractor.SIGNATURES.contains(featureName)) {
            return FeatureType.SIGNATURE;
        } else if (DetermineCopyNumber.isAmplification(featureName, biomarkerType)) {
            return FeatureType.AMPLIFICATION;
        } else if (DetermineCopyNumber.isDeletion(featureName, biomarkerType) && !featureName.toLowerCase().contains("exon")) {
            return FeatureType.DELETION;
        } else if (DetermineFusion.isFusion(featureName, biomarkerType, provenanceRule, proteinAnnotation) && !featureName.contains(
                "p61BRAF")) {
            return FeatureType.FUSION_PAIR;
        } else if (DetermineFusion.isFusionPromiscuous(featureName, biomarkerType, provenanceRule, proteinAnnotation)) {
            return FeatureType.FUSION_PROMISCUOUS;
        } else if (FeatureTypeExtractor.GENE_EXON.contains(event)) {
            if (featureName.toLowerCase().contains("deletion") || featureName.toLowerCase().contains("insertion")
                    || featureName.toLowerCase().contains("proximal") || featureName.toLowerCase().contains("mutation")
                    || featureName.toLowerCase().contains("splice site insertion") || featureName.toLowerCase().contains("frameshift")) {
                return FeatureType.GENE_RANGE_EXON;
            }
        } else if (proteinAnnotation.length() > 1 && proteinAnnotation.substring(proteinAnnotation.length() - 1).equals("X")) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (proteinAnnotation.length() >= 1 && isValidSingleCodonRange(proteinAnnotation)) {
            return FeatureType.GENE_RANGE_CODON;
        } else if (!DetermineHotspot.isHotspot(proteinAnnotation)) {
            String eventDescription = featureDescription.split(" ", 2)[1].trim();
            if (FeatureTypeExtractor.DETAILED_GENE_LEVEL_INFO_WITHOUT_TSG_ONCO.contains(eventDescription)
                    || FeatureTypeExtractor.DETAILED_GENE_LEVEL_INFO_WITH_TSG.contains(eventDescription)
                    || FeatureTypeExtractor.DETAILED_GENE_LEVEL_INFO_WITH_ONCO.contains(eventDescription) || FeatureTypeExtractor.GENE_LEVEL
                    .contains(provenanceRule)) {
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
