package com.hartwig.hmftools.common.serve;

import java.util.Comparator;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum Knowledgebase {
    DOCM("docm"),
    HARTWIG_COHORT("hartwig_cohort"),
    HARTWIG_CURATED("hartwig_curated"),
    ICLUSION("iclusion"),
    VICC_CGI("vicc_cgi"),
    VICC_CIVIC("vicc_civic"),
    VICC_JAX("vicc_jax"),
    VICC_ONCOKB("vicc_oncokb");

    @NotNull
    private final String display;

    private static final Logger LOGGER = LogManager.getLogger(Knowledgebase.class);

    Knowledgebase(@NotNull final String display) {
        this.display = display;
    }

    @NotNull
    public String display() {
        return display;
    }

    @NotNull
    public static Set<Knowledgebase> fromCommaSeparatedSourceString(@NotNull String sources) {
        Set<Knowledgebase> consolidated = Sets.newHashSet();

        for (String source : sources.split(",")) {
            Knowledgebase knowledgebase = lookupKnowledgebase(source);
            if (knowledgebase != null) {
                consolidated.add(knowledgebase);
            } else {
                LOGGER.warn("Could not resolve knowledgebase with display '{}'", source);
            }
        }

        return consolidated;
    }

    @Nullable
    public static Knowledgebase lookupKnowledgebase(@NotNull String display) {
        for (Knowledgebase knowledgebase : Knowledgebase.values()) {
            if (knowledgebase.display().equals(display)) {
                return knowledgebase;
            }
        }
        return null;
    }

    @NotNull
    public static String toCommaSeparatedSourceString(@NotNull Set<Knowledgebase> sources) {
        Set<Knowledgebase> sorted = Sets.newTreeSet(Comparator.naturalOrder());
        sorted.addAll(sources);

        StringJoiner joiner = new StringJoiner(",");
        for (Knowledgebase source : sorted) {
            joiner.add(source.toString().toLowerCase());
        }
        return joiner.toString();
    }
}
