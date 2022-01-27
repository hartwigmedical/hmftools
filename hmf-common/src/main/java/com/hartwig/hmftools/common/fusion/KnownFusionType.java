package com.hartwig.hmftools.common.fusion;

import org.jetbrains.annotations.NotNull;

public enum KnownFusionType
{
    NONE("None"),
    KNOWN_PAIR("Known pair"),
    PROMISCUOUS_5("5' Promiscuous"),
    PROMISCUOUS_3("3' Promiscuous"),
    IG_KNOWN_PAIR("IG known pair"),
    IG_PROMISCUOUS("IG promiscuous"),
    EXON_DEL_DUP("Exon del dup"),
    PROMISCUOUS_BOTH("5' and 3' Promiscuous");

    private final String mDisplay;

    KnownFusionType(@NotNull final String display) {
        mDisplay = display;
    }

    public String display() {
        return mDisplay;
    }
}
