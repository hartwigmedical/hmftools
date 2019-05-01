package com.hartwig.hmftools.patientreporter.variants.germline;

import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GermlineReportingModel {

    private static final Logger LOGGER = LogManager.getLogger(GermlineReportingModel.class);

    @NotNull
    private final Map<String, Boolean> germlineGenesAndAndNotificationMap;

    GermlineReportingModel(@NotNull final Map<String, Boolean> germlineGenesAndAndNotificationMap) {
        this.germlineGenesAndAndNotificationMap = germlineGenesAndAndNotificationMap;
    }

    @NotNull
    public Set<String> reportableGermlineGenes() {
        return germlineGenesAndAndNotificationMap.keySet();
    }

    public boolean genesToNotifyClinicalGeneticist(@NotNull String germlineGene) {
        Boolean notify = germlineGenesAndAndNotificationMap.get(germlineGene);
        if (notify == null) {
            LOGGER.warn("Requested notification status for a gene that is not amongst set of reportable germline genes: " + germlineGene);
        }

        return notify != null ? notify : false;
    }
}
