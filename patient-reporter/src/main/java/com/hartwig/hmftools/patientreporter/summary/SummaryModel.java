package com.hartwig.hmftools.patientreporter.summary;

import java.util.Map;
import java.util.Set;

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
    public Set<String> sizeSummarySamples() {
        return sampleToSummaryMap.keySet();
    }

    public boolean sampleIdPresentInSummaryFile(@NotNull String sampleId) {
        return sampleToSummaryMap.keySet().contains(sampleId);
    }

    @NotNull
    public String extractSummarySampleId(@NotNull String sampleId) {
        boolean sampleInFile = sampleIdPresentInSummaryFile(sampleId);
        if (!sampleInFile) {
            LOGGER.info("Could not find a summary for: " + sampleId);
        }
        return sampleInFile ? sampleToSummaryMap.get(sampleId) : Strings.EMPTY;
    }
}
