package com.hartwig.hmftools.patientreporter.germline;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GermlineReportingModel {

    private static final Logger LOGGER = LogManager.getLogger(GermlineReportingModel.class);

    @NotNull
    private final List<GermlineReportingEntry> entries;

    public GermlineReportingModel(@NotNull final List<GermlineReportingEntry> entries) {
        this.entries = entries;
    }

    @NotNull
    @VisibleForTesting
    List<GermlineReportingEntry> entries() {
        return entries;
    }

    @Nullable
    public GermlineReportingEntry entryForGene(@NotNull String gene) {
        for (GermlineReportingEntry entry : entries) {
            if (entry.gene().equals(gene)) {
                return entry;
            }
        }
        return null;
    }

    public boolean notifyAboutGene(@NotNull String gene, @NotNull LimsGermlineReportingLevel reportingLevel) {
        GermlineReportingEntry entry = entryForGene(gene);
        if (entry == null) {
            LOGGER.warn("Requested notification status for a gene that is not amongst set of reportable germline genes: {}", gene);
            return false;
        }

        return reportingLevel == LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION && entry.notifyClinicalGeneticist();
    }
}
