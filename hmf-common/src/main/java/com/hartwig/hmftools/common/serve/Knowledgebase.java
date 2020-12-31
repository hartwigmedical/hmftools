package com.hartwig.hmftools.common.serve;

import java.util.Comparator;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum Knowledgebase {
    DOCM("docm"),
    HARTWIG_COHORT("hartwig_cohort"),
    HARTWIG_CURATED("hartwig_curated"),
    ICLUSION("iclusion"),
    VICC_CGI("vicc_cgi"),
    VICC_CIVIC("vicc_civic"),
    VICC_JAX("vicc_jax"),
    VICC_ONCOKB("vicc_oncokb"),
    UNDEFINED("undefined");

    private static final Logger LOGGER = LogManager.getLogger(Knowledgebase.class);

    @NotNull
    private final String knowledgebase;

    Knowledgebase(@NotNull final String knowledgebase) {
        this.knowledgebase = knowledgebase;
    }

    @NotNull
    public String knowledgebase() {
        return knowledgebase;
    }

    @NotNull
    public static Set<Knowledgebase> extractKnowledgebase(@NotNull String knowledgebase) {
        Set<Knowledgebase> consolidated = Sets.newHashSet();

        String [] multipleSources = knowledgebase.split(",");
        for (String sources: multipleSources) {

            for (Knowledgebase knowledgebaseSource : Knowledgebase.values()) {
                if (knowledgebaseSource.knowledgebase().equals(sources)) {
                    consolidated.add(knowledgebaseSource);
                }
            }
        }

        return consolidated;

    }

    @NotNull
    public static String commaSeparatedSourceString(@NotNull Set<Knowledgebase> sources) {
        Set<Knowledgebase> sorted = Sets.newTreeSet(Comparator.naturalOrder());
        sorted.addAll(sources);

        StringJoiner joiner = new StringJoiner(",");
        for (Knowledgebase source : sorted) {
            joiner.add(source.toString().toLowerCase());
        }
        return joiner.toString();
    }
}
