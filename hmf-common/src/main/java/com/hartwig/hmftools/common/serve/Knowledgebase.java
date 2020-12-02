package com.hartwig.hmftools.common.serve;

import java.util.Comparator;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public enum Knowledgebase {
    DOCM,
    HARTWIG_COHORT,
    HARTWIG_CURATED,
    ICLUSION,
    VICC_CGI,
    VICC_CIVIC,
    VICC_JAX,
    VICC_ONCOKB;

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
