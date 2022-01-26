package com.hartwig.hmftools.common.fusion;

import org.jetbrains.annotations.NotNull;

public enum KnownFusionType
{
    NONE("none"),
    KNOWN_PAIR("known pair"),
    PROMISCUOUS_5("5' promiscuous"),
    PROMISCUOUS_3("3' promiscuous"),
    IG_KNOWN_PAIR("ig known pair"),
    IG_PROMISCUOUS("ig promiscuous"),
    EXON_DEL_DUP("exon del dup"),
    PROMISCUOUS_BOTH("5' en 3' promiscuous");

    private final String mDisplay;

    KnownFusionType(@NotNull final String display) {
        mDisplay = display;
    }

    public String display() {
        return mDisplay;
    }
}
