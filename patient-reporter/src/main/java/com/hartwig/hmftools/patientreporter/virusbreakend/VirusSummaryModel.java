package com.hartwig.hmftools.patientreporter.virusbreakend;

import java.util.Collection;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class VirusSummaryModel {

    private static final Logger LOGGER = LogManager.getLogger(VirusSummaryModel.class);

    @NotNull
    private final Map<Integer, String> virusSummaryMap;

    public VirusSummaryModel(@NotNull final Map<Integer, String> virusSummaryMap) {
        this.virusSummaryMap = virusSummaryMap;
    }

    @NotNull
    public String findVirusSummary(int id) {
        boolean mappedVirusSummary = mapIdtoVirusName(id);

        if (!mappedVirusSummary) {
            LOGGER.warn("Could not match id to virusName");
        }
        return mappedVirusSummary ? virusSummaryMap.get(id) : Strings.EMPTY;
    }

    @VisibleForTesting
    public boolean mapIdtoVirusName(int id) {
        return virusSummaryMap.containsKey(id);
    }

    @VisibleForTesting
    public Collection<String> virussen() {
        return virusSummaryMap.values();
    }
}
