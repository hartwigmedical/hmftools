package com.hartwig.hmftools.knowledgebasegenerator.fusion;

import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class FusionExtractor {
    private static final Logger LOGGER = LogManager.getLogger(FusionExtractor.class);

    @NotNull
    public static KnownFusions determineKnownFusions(@NotNull Source source, @NotNull String typeEvent,
            @NotNull String gene) {
        return ImmutableKnownFusions.builder().gene(gene).eventType(typeEvent).source(source.toString()).sourceLink(source.toString()).build();

    }

    @NotNull
    public static KnownFusions determineKnownFusionsPairs(@NotNull Source source, @NotNull String typeEvent,
            @NotNull String gene) {
        return ImmutableKnownFusions.builder().gene(gene).eventType(typeEvent).source(source.toString()).sourceLink(source.toString()).build();

    }

}
