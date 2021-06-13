package com.hartwig.hmftools.common.serve;

import java.util.Comparator;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public enum Knowledgebase {
    CKB(RefGenomeVersion.V38, "ckb", "CKB", EvidenceLevel.C, EvidenceLevel.B),
    DOCM(RefGenomeVersion.V37, "docm", "DoCM", EvidenceLevel.A, EvidenceLevel.A),
    HARTWIG_COHORT(RefGenomeVersion.V37, "hartwig_cohort", "HMF Cohort", EvidenceLevel.A, EvidenceLevel.A),
    HARTWIG_CURATED(RefGenomeVersion.V37, "hartwig_curated", "HMF Curated", EvidenceLevel.A, EvidenceLevel.A),
    ICLUSION(RefGenomeVersion.V37, "iclusion", "iClusion", EvidenceLevel.B, EvidenceLevel.B),
    VICC_CGI(RefGenomeVersion.V37, "vicc_cgi", "CGI", EvidenceLevel.B, EvidenceLevel.B),
    VICC_CIVIC(RefGenomeVersion.V37, "vicc_civic", "CIViC", EvidenceLevel.B, EvidenceLevel.B),
    VICC_JAX(RefGenomeVersion.V37, "vicc_jax", "CKB Core", EvidenceLevel.B, EvidenceLevel.B),
    VICC_ONCOKB(RefGenomeVersion.V37, "vicc_oncokb", "OncoKB", EvidenceLevel.B, EvidenceLevel.B);

    private static final Logger LOGGER = LogManager.getLogger(Knowledgebase.class);

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final String technicalDisplay;
    @NotNull
    private final String reportDisplay;
    @NotNull
    private final EvidenceLevel maxCertainEvidenceReportingLevel;
    @NotNull
    private final EvidenceLevel maxPredictedEvidenceReportingLevel;

    Knowledgebase(@NotNull final RefGenomeVersion refGenomeVersion, @NotNull final String technicalDisplay,
            @NotNull final String reportDisplay, @NotNull final EvidenceLevel maxCertainEvidenceReportingLevel,
            @NotNull final EvidenceLevel maxPredictedEvidenceReportingLevel) {
        this.refGenomeVersion = refGenomeVersion;
        this.technicalDisplay = technicalDisplay;
        this.reportDisplay = reportDisplay;
        this.maxCertainEvidenceReportingLevel = maxCertainEvidenceReportingLevel;
        this.maxPredictedEvidenceReportingLevel = maxPredictedEvidenceReportingLevel;
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
    public EvidenceLevel maxCertainEvidenceReportingLevel() {
        return maxCertainEvidenceReportingLevel;
    }

    @NotNull
    public EvidenceLevel maxPredictedEvidenceReportingLevel() {
        return maxPredictedEvidenceReportingLevel;
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
