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

public enum Knowledgebase {
    CKB(RefGenomeVersion.V38, "CKB", EvidenceLevel.C, EvidenceLevel.B),
    DOCM(RefGenomeVersion.V37, "DoCM", EvidenceLevel.A, EvidenceLevel.A),
    HARTWIG_COHORT(RefGenomeVersion.V37, "HMF Cohort", EvidenceLevel.A, EvidenceLevel.A),
    HARTWIG_CURATED(RefGenomeVersion.V37, "HMF Curated", EvidenceLevel.A, EvidenceLevel.A),
    ICLUSION(RefGenomeVersion.V37, "iClusion", EvidenceLevel.B, EvidenceLevel.B),
    ACTIN(RefGenomeVersion.V37, "ACTIN", EvidenceLevel.B, EvidenceLevel.B),
    VICC_CGI(RefGenomeVersion.V37, "CGI", EvidenceLevel.B, EvidenceLevel.B),
    VICC_CIVIC(RefGenomeVersion.V37, "CIViC", EvidenceLevel.B, EvidenceLevel.B),
    VICC_JAX(RefGenomeVersion.V37, "CKB Core", EvidenceLevel.B, EvidenceLevel.B),
    VICC_ONCOKB(RefGenomeVersion.V37, "OncoKB", EvidenceLevel.B, EvidenceLevel.B),
    UNKNOWN(RefGenomeVersion.V37, "Unknown", EvidenceLevel.D, EvidenceLevel.D);

    private static final Logger LOGGER = LogManager.getLogger(Knowledgebase.class);

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final String display;
    @NotNull
    private final EvidenceLevel maxCertainEvidenceReportingLevel;
    @NotNull
    private final EvidenceLevel maxPredictedEvidenceReportingLevel;

    Knowledgebase(@NotNull final RefGenomeVersion refGenomeVersion, @NotNull final String display,
            @NotNull final EvidenceLevel maxCertainEvidenceReportingLevel,
            @NotNull final EvidenceLevel maxPredictedEvidenceReportingLevel) {
        this.refGenomeVersion = refGenomeVersion;
        this.display = display;
        this.maxCertainEvidenceReportingLevel = maxCertainEvidenceReportingLevel;
        this.maxPredictedEvidenceReportingLevel = maxPredictedEvidenceReportingLevel;
    }

    @NotNull
    public RefGenomeVersion refGenomeVersion() {
        return refGenomeVersion;
    }

    @NotNull
    public String display() {
        return display;
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
    public static Set<Knowledgebase> fromCommaSeparatedSourceString(@NotNull String sources) {
        Set<Knowledgebase> consolidated = Sets.newHashSet();

        for (String source : sources.split(",")) {
            Knowledgebase knowledgebase = lookupKnowledgebase(source);
            if (knowledgebase != Knowledgebase.UNKNOWN) {
                consolidated.add(knowledgebase);
            } else {
                LOGGER.warn("Could not resolve knowledgebase from source '{}'", source);
            }
        }

        return consolidated;
    }

    @NotNull
    public static Knowledgebase lookupKnowledgebase(@NotNull String knowledgebaseToFind) {
        for (Knowledgebase knowledgebase : Knowledgebase.values()) {
            if (knowledgebase.toString().equals(knowledgebaseToFind)) {
                return knowledgebase;
            }
        }
        return Knowledgebase.UNKNOWN;
    }

    @NotNull
    public static String toCommaSeparatedSourceString(@NotNull Set<Knowledgebase> sources) {
        Set<Knowledgebase> sorted = Sets.newTreeSet(Comparator.naturalOrder());
        sorted.addAll(sources);

        StringJoiner joiner = new StringJoiner(",");
        for (Knowledgebase source : sorted) {
            joiner.add(source.toString());
        }
        return joiner.toString();
    }
}