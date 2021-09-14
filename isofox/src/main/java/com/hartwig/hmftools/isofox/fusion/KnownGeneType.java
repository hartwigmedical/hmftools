package com.hartwig.hmftools.isofox.fusion;

public enum KnownGeneType
{
    KNOWN_PAIR,
    KNOWN_PROM3,
    KNOWN_OTHER,
    PROM5_KNOWN,
    PROM5_PROM3,
    PROM5_OTHER,
    OTHER_PROM3,
    OTHER;

    public static boolean hasKnownPairGene(final KnownGeneType type)
    {
        return type == KNOWN_PAIR || type == KNOWN_PROM3 || type == KNOWN_OTHER || type == PROM5_KNOWN;
    }
}
