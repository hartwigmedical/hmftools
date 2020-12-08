package com.hartwig.hmftools.patientreporter.summary;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.lims.LimsCohort;

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
    public String findSummaryForSample(@NotNull String sample, @NotNull LimsCohort cohort) {
        boolean sampleHasSummary = samplePresentInSummaries(sample);

        if (!sampleHasSummary && cohort == LimsCohort.WIDE) {
            LOGGER.warn("Could not find a summary for WIDE sample: {}", sample);
        } else if (cohort == LimsCohort.CORE || cohort == LimsCohort.CORELR11 || cohort == LimsCohort.CORESC11 && !sampleHasSummary) {
            LOGGER.warn("Could not find a summary for CORE sample: {}", sample);
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
