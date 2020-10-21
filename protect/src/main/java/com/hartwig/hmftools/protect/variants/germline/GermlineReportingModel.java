package com.hartwig.hmftools.protect.variants.germline;

import java.util.Collections;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GermlineReportingModel {

    private static final Logger LOGGER = LogManager.getLogger(GermlineReportingModel.class);

    @NotNull
    private final Set<String> reportableGenes = Sets.newHashSet();
    private final Set<String> notifiableGenes = Sets.newHashSet();;

    public GermlineReportingModel(@NotNull final Map<String, Boolean> germlineGenesAndNotificationMap) {
        for (Map.Entry<String, Boolean> entry : germlineGenesAndNotificationMap.entrySet()) {
            reportableGenes.add(entry.getKey());
            if (entry.getValue()) {
                notifiableGenes.add(entry.getKey());
            }
        }
    }

    @NotNull
    public Set<String> reportableGermlineGenes() {
        return reportableGenes;
    }

    @NotNull
    Set<String> notifiableGenes(@NotNull LimsGermlineReportingLevel reportingLevel) {
        if (reportingLevel.equals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION)) {
            return notifiableGenes;
        }

        return Collections.emptySet();
    }

    public boolean notifyAboutGene(@NotNull LimsGermlineReportingLevel reportingLevel, @NotNull String germlineGene) {
        boolean reportableContains = reportableGenes.contains(germlineGene);
        if (!reportableContains) {
            LOGGER.warn("Requested notification status for a gene that is not amongst set of reportable germline genes: {}", germlineGene);
        }

        return reportingLevel.equals(LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION) &&  notifiableGenes.contains(germlineGene);
    }
}
