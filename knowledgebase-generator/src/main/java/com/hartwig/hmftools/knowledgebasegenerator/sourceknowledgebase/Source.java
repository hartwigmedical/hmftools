package com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase;

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
            return UNKNOWN;
        }
    }

}
