package com.hartwig.hmftools.vicc.datamodel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum ViccSource {
    BRCA("brca"),
    CGI("cgi"),
    CIVIC("civic"),
    JAX("jax"),
    JAX_TRIALS("jax_trials"),
    MOLECULAR_MATCH("molecularmatch"),
    MOLECULAR_MATCH_TRIALS("molecularmatch_trials"),
    ONCOKB("oncokb"),
    PMKB("pmkb"),
    SAGE("sage"),
    UNKNOWN("unknown");

    private static final Logger LOGGER = LogManager.getLogger(ViccSource.class);

    @NotNull
    private final String display;

    ViccSource(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }

    @NotNull
    public static ViccSource fromViccKnowledgebaseString(@NotNull String sourceString) {
        for (ViccSource source : ViccSource.values()) {
            if (source.display().equals(sourceString)) {
                return source;
            }
        }

        LOGGER.warn("Unknown source in knowledgebase: {}!", sourceString);
        return ViccSource.UNKNOWN;
    }
}
