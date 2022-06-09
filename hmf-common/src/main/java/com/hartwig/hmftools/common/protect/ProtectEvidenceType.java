package com.hartwig.hmftools.common.protect;

import org.jetbrains.annotations.NotNull;

public enum ProtectEvidenceType {
    VIRAL_PRESENCE("Viral"),
    SIGNATURE("Signature"),
    ACTIVATION("Activation"),
    INACTIVATION("Inactivation"),
    AMPLIFICATION("Amplification"),
    OVER_EXPRESSION("Over expression"),
    DELETION("Deletion"),
    UNDER_EXPRESSION("Under expression"),
    PROMISCUOUS_FUSION("Promiscuous fusion"),
    FUSION_PAIR("Fusion pair"),
    HOTSPOT_MUTATION("Hotspot"),
    CODON_MUTATION("Codon"),
    EXON_MUTATION("Exon"),
    ANY_MUTATION("Any mutation"),
    WILD_TYPE("Wild-type"),
    HLA("hla");

    private final String mDisplay;

    ProtectEvidenceType(@NotNull final String display) {
        mDisplay = display;
    }

    public String display() {
        return mDisplay;
    }
}