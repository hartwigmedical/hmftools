package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LimsModel {

    private static final Logger LOGGER = LogManager.getLogger(LimsModel.class);

    @NotNull
    private final Map<String, LimsData> dataPerSample;

    LimsModel(@NotNull final Map<String, LimsData> dataPerSample) {
        this.dataPerSample = dataPerSample;
    }

    @NotNull
    @VisibleForTesting
    Map<String, LimsData> data() {
        return dataPerSample;
    }

    @Nullable
    public Double findTumorPercentageForSample(@NotNull final String sample) {
        final LimsData dataForSample = dataPerSample.get(sample);
        if (dataForSample == null) {
            LOGGER.warn("Could not find LIMS data for " + sample);
            return null;
        }
        if (dataForSample instanceof LimsBloodData) {
            LOGGER.warn("Sample " + sample + " is a blood sample and has no tumor percentage.");
            return null;
        }
        if (dataForSample instanceof LimsTumorData) {
            final LimsTumorData tumorSample = (LimsTumorData) dataForSample;
            return tumorSample.tumorPercentage();
        }
        return null;
    }

    @Nullable
    public LocalDate findArrivalDateForSample(@NotNull final String sample) {
        final LimsData dataForSample = dataPerSample.get(sample);
        if (dataForSample == null) {
            LOGGER.warn(" Could not find LIMS data for " + sample);
            return null;
        }
        return dataForSample.arrivalDate();
    }

    @Nullable
    public LocalDate findSamplingDateForSample(@NotNull final String sample) {
        final LimsData dataForSample = dataPerSample.get(sample);
        if (dataForSample == null) {
            LOGGER.warn(" Could not find LIMS data for " + sample);
            return null;
        }
        return dataForSample.samplingDate();
    }

    @Nullable
    public LimsData findDataPerSample(@NotNull final String sample) {
        return dataPerSample.get(sample);
    }
}
