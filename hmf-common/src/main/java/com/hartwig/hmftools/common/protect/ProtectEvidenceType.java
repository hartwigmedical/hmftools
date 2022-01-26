package com.hartwig.hmftools.common.protect;

import org.jetbrains.annotations.NotNull;

public enum ProtectEvidenceType {
    VIRAL_PRESENCE("viral presence"),
    SIGNATURE("signature"),
    ACTIVATION("activation"),
    INACTIVATION("inactivation"),
    AMPLIFICATION("amplification"),
    DELETION("deletion"),
    PROMISCUOUS_FUSION("promiscuous fusion"),
    FUSION_PAIR("fusion pair"),
    HOTSPOT_MUTATION("hotspot"),
    CODON_MUTATION("codon"),
    EXON_MUTATION("exon"),
    ANY_MUTATION("any mutation");

    private final String mDisplay;

    ProtectEvidenceType(@NotNull final String display) {
        mDisplay = display;
    }

    public String display() {
        return mDisplay;
    }
}