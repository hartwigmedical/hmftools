package com.hartwig.hmftools.patientreporter.summary;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.lims.LimsCohort;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SummaryModel {

    private static final Logger LOGGER = LogManager.getLogger(SummaryModel.class);

    @NotNull
    private final Map<String, String> sampleToSummaryMap;

    SummaryModel(@NotNull final Map<String, String> sampleToSummaryMap) {
        this.sampleToSummaryMap = sampleToSummaryMap;
    }

    @NotNull
    public String findSummaryForSample(@NotNull String sample, @NotNull LimsCohortConfigData cohort) {
        boolean sampleHasSummary = samplePresentInSummaries(sample);

        if (!sampleHasSummary && cohort.reportConclusion()) {
            LOGGER.warn("Could not find a summary for WIDE/CORE sample: {}", sample);
        }
        return sampleHasSummary ? sampleToSummaryMap.get(sample) : Strings.EMPTY;
    }

    @VisibleForTesting
    int summaryCount() {
        return sampleToSummaryMap.keySet().size();
    }

    @VisibleForTesting
    boolean samplePresentInSummaries(@NotNull String sample) {
        return sampleToSummaryMap.containsKey(sample);
    }
}
