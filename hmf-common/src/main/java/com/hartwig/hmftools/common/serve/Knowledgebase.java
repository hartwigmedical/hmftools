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
    DOCM("docm", "DoCM"),
    HARTWIG_COHORT("hartwig_cohort", "HMF Cohort"),
    HARTWIG_CURATED("hartwig_curated", "HMF Curated"),
    ICLUSION("iclusion", "iClusion"),
    VICC_CGI("vicc_cgi", "CGI"),
    VICC_CIVIC("vicc_civic", "CIViC"),
    VICC_JAX("vicc_jax", "CKB"),
    VICC_ONCOKB("vicc_oncokb", "OncoKB");

    @NotNull
    private final String technicalDisplay;
    @NotNull
    private final String reportDisplay;

    private static final Logger LOGGER = LogManager.getLogger(Knowledgebase.class);

    Knowledgebase(@NotNull final String technicalDisplay, @NotNull final String reportDisplay) {
        this.technicalDisplay = technicalDisplay;
        this.reportDisplay = reportDisplay;
    }

    @NotNull
    public String technicalDisplay() {
        return technicalDisplay;
    }

    @NotNull
    public String reportDisplay() {
        return reportDisplay;
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
    private static Knowledgebase lookupKnowledgebase(@NotNull String technicalDisplay) {
        for (Knowledgebase knowledgebase : Knowledgebase.values()) {
            if (knowledgebase.technicalDisplay().equals(technicalDisplay)) {
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
            joiner.add(source.technicalDisplay());
        }
        return joiner.toString();
    }
}
