package com.hartwig.hmftools.common.rna;

public enum KnownFusionType
{
    KNOWN_PAIR,
    KNOWN_PROM3,
    KNOWN_OTHER,
    PROM5_KNOWN,
    PROM5_PROM3,
    PROM5_OTHER,
    OTHER_PROM3,
    OTHER;

    public static boolean hasKnownPairGene(final KnownFusionType type)
    {
        return type == KNOWN_PAIR || type == KNOWN_PROM3 || type == KNOWN_OTHER || type == PROM5_KNOWN;
    }

    public static boolean hasPromiscousGene(final KnownFusionType type)
    {
        return type == PROM5_PROM3 || type == PROM5_KNOWN || type == PROM5_OTHER
            ||  type == KNOWN_PROM3 || type == OTHER_PROM3;
    }
}
