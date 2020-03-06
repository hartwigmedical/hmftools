package com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase;

import com.hartwig.hmftools.knowledgebasegenerator.eventtype.EventTypeAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum Source {
    ONCOKB,
    CGI,
    CIVIC,
    JAX,
    JAX_TRIALS,
    BRCA,
    SAGE,
    PMKB,
    MOLECULARMATCH,
    MOLECULARMATCH_TRIALS,
    UNKNOWN;

    private static final Logger LOGGER = LogManager.getLogger(EventTypeAnalyzer.class);

    @NotNull
    public static Source sourceFromKnowledgebase(@NotNull String source) {
        if (source.equals("oncokb")) {
            return ONCOKB;
        } else if (source.equals("cgi")) {
            return CGI;
        } else if (source.equals("civic")) {
            return CIVIC;
        } else if (source.equals("jax")) {
            return JAX;
        } else if (source.equals("jax_trials")) {
            return JAX_TRIALS;
        } else if (source.equals("brca")) {
            return BRCA;
        } else if (source.equals("sage")) {
            return SAGE;
        } else if (source.equals("pmkb")) {
            return PMKB;
        } else if (source.equals("molecularmatch")) {
            return MOLECULARMATCH;
        } else if (source.equals("molecularmatch_trials")) {
            return MOLECULARMATCH_TRIALS;
        } else {
            LOGGER.warn("Unknown source in knowledgebase!");
            return UNKNOWN;
        }
    }

}
