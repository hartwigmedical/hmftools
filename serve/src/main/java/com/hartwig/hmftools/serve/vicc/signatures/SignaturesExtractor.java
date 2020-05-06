package com.hartwig.hmftools.serve.vicc.signatures;

import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SignaturesExtractor {

    private static final Logger LOGGER = LogManager.getLogger(SignaturesExtractor.class);

    @NotNull
    public static Signatures determineSignatures(@NotNull ViccSource source, @NotNull String typeEvent) {
        return ImmutableSignatures.builder().eventType(typeEvent).source(source.displayString()).sourceLink(source.toString()).build();
    }
}
