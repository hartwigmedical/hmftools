package com.hartwig.hmftools.vicc.util;

import java.util.Set;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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

    public static final Set<String> SIGNATURES = Sets.newHashSet("Microsatellite Instability-High");

    @NotNull
    public static EventAnnotation toEventAnnotation(@NotNull String featureName, @Nullable String biomarkerType) {
        String feature = featureName;
        if (feature.contains(" ") && !feature.equals("Copy Number Loss")) {
            feature = feature.split(" ", 2)[1];
        }

        if (EventAnnotationExtractor.SIGNATURES.contains(feature)) {
            return EventAnnotation.SIGNATURE;
        } else if (EventAnnotationExtractor.AMPLIFICATIONS.contains(feature) ||
                EventAnnotationExtractor.AMPLIFICATIONS.contains(biomarkerType)) {
            return EventAnnotation.AMPLIFICATION;
        } else if (EventAnnotationExtractor.DELETIONS.contains(feature) ||
                EventAnnotationExtractor.DELETIONS.contains(biomarkerType)) {
            return EventAnnotation.DELETION;
        } else {
//            LOGGER.warn("No event annotation extracted from event!");
            return EventAnnotation.UNKNOWN;
        }
    }
}
