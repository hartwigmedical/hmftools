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
    CKB(RefGenomeVersion.V38, "ckb", "CKB"),
    DOCM(RefGenomeVersion.V37, "docm", "DoCM"),
    HARTWIG_COHORT(RefGenomeVersion.V37, "hartwig_cohort", "HMF Cohort"),
    HARTWIG_CURATED(RefGenomeVersion.V37, "hartwig_curated", "HMF Curated"),
    ICLUSION(RefGenomeVersion.V37, "iclusion", "iClusion"),
    VICC_CGI(RefGenomeVersion.V37, "vicc_cgi", "CGI"),
    VICC_CIVIC(RefGenomeVersion.V37, "vicc_civic", "CIViC"),
    VICC_JAX(RefGenomeVersion.V37, "vicc_jax", "CKB Core"),
    VICC_ONCOKB(RefGenomeVersion.V37, "vicc_oncokb", "OncoKB");

    private static final Logger LOGGER = LogManager.getLogger(Knowledgebase.class);

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final String technicalDisplay;
    @NotNull
    private final String reportDisplay;

    Knowledgebase(@NotNull final RefGenomeVersion refGenomeVersion, @NotNull final String technicalDisplay,
            @NotNull final String reportDisplay) {
        this.refGenomeVersion = refGenomeVersion;
        this.technicalDisplay = technicalDisplay;
        this.reportDisplay = reportDisplay;
    }

    @NotNull
    public RefGenomeVersion refGenomeVersion() {
        return refGenomeVersion;
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
    public static Set<Knowledgebase> fromCommaSeparatedTechnicalDisplayString(@NotNull String knowledgebases) {
        Set<Knowledgebase> consolidated = Sets.newHashSet();

        for (String technicalDisplay : knowledgebases.split(",")) {
            Knowledgebase knowledgebase = lookupKnowledgebase(technicalDisplay);
            if (knowledgebase != null) {
                consolidated.add(knowledgebase);
            } else {
                LOGGER.warn("Could not resolve knowledgebase with display '{}'", knowledgebase);
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
