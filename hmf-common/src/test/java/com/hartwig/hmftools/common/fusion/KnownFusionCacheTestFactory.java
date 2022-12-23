package com.hartwig.hmftools.common.fusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class KnownFusionCacheTestFactory {

    private KnownFusionCacheTestFactory()
    {
    }

    @NotNull
    public static KnownFusionData createKnownPair(@NotNull String fiveGene, @NotNull String threeGene)
    {
        return create(KnownFusionType.KNOWN_PAIR, fiveGene, threeGene);
    }

    @NotNull
    public static KnownFusionData createPromiscuousFive(@NotNull String gene)
    {
        return create(KnownFusionType.PROMISCUOUS_5, gene, Strings.EMPTY);
    }

    @NotNull
    public static KnownFusionData createPromiscuousThree(@NotNull String gene)
    {
        return create(KnownFusionType.PROMISCUOUS_3, Strings.EMPTY, gene);
    }

    @NotNull
    public static KnownFusionData createExonDelDup(@NotNull String gene)
    {
        return create(KnownFusionType.EXON_DEL_DUP, gene, gene);
    }

    @NotNull
    public static KnownFusionData create(@NotNull KnownFusionType type, @NotNull String fiveGene, @NotNull String threeGene)
    {
        return new KnownFusionData(type, fiveGene, threeGene, Strings.EMPTY, Strings.EMPTY);
    }
}
