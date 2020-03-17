package com.hartwig.hmftools.knowledgebasegenerator.signatures;

import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {

    private static final Logger LOGGER = LogManager.getLogger(SignaturesExtractor.class);

    @NotNull
    public static Signatures determineSignatures(@NotNull Source source, @NotNull String typeEvent) {
        return ImmutableSignatures.builder().eventType(typeEvent).source(source.toString()).sourceLink(source.toString()).build();

    }

}
