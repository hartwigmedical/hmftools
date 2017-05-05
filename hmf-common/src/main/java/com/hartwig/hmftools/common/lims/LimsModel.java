package com.hartwig.hmftools.common.lims;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LimsModel {

    private static final Logger LOGGER = LogManager.getLogger(LimsModel.class);

    @NotNull
    private final Map<String, LimsBiopsyData> dataPerSample;

    LimsModel(@NotNull final Map<String, LimsBiopsyData> dataPerSample) {
        this.dataPerSample = dataPerSample;
    }

    @NotNull
    @VisibleForTesting
    Map<String, LimsBiopsyData> data() {
        return dataPerSample;
    }

    public double findTumorPercentageForSample(@NotNull final String sample) {
        final LimsBiopsyData dataForSample = dataPerSample.get(sample);
        if (dataForSample == null) {
            LOGGER.warn(" Could not find LIMS data for " + sample);
            return Double.NaN;
        }
        return dataForSample.tumorPercentage();
    }

    public String findArrivalDateForSample(@NotNull final String sample) {
        final LimsBiopsyData dataForSample = dataPerSample.get(sample);
        if (dataForSample == null) {
            LOGGER.warn(" Could not find LIMS data for " + sample);
            return null;
        }
        return dataForSample.arrivalDate();
    }

    public String findSamplingDateForSample(@NotNull final String sample) {
        final LimsBiopsyData dataForSample = dataPerSample.get(sample);
        if (dataForSample == null) {
            LOGGER.warn(" Could not find LIMS data for " + sample);
            return null;
        }
        return dataForSample.samplingDate();
    }
}
