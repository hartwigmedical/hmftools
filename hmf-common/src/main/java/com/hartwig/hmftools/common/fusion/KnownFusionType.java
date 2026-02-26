package com.hartwig.hmftools.common.fusion;

public enum KnownFusionType
{
    NONE,
    KNOWN_PAIR,
    PROMISCUOUS_5,
    PROMISCUOUS_3,
    ENHANCER_KNOWN_PAIR,
    ENHANCER_PROMISCUOUS,
    EXON_DEL_DUP,
    PROMISCUOUS_ENHANCER_TARGET;

    public static final String PROMISCUOUS_BOTH = "PROMISCUOUS_BOTH";

    // pre v3.0
    public static final String IG_KNOWN_PAIR_STR = "IG_KNOWN_PAIR";
    public static final String IG_PROMISCUOUS_STR = "IG_PROMISCUOUS";

    public static KnownFusionType parse(final String knownFusionTypeStr)
    {
        if(knownFusionTypeStr.equals(IG_KNOWN_PAIR_STR))
            return ENHANCER_KNOWN_PAIR;

        if(knownFusionTypeStr.equals(IG_PROMISCUOUS_STR))
            return ENHANCER_PROMISCUOUS;

        return KnownFusionType.valueOf(knownFusionTypeStr);
    }

}
