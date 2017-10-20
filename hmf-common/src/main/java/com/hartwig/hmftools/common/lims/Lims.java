package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Lims {

    private static final Logger LOGGER = LogManager.getLogger(LimsJsonModel.class);

    @VisibleForTesting
    static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final Map<String, LimsJsonData> dataPerSample;
    @NotNull
    private final Map<String, LocalDate> preHmfArrivalDates;

    Lims(@NotNull final Map<String, LimsJsonData> dataPerSample, @NotNull final Map<String, LocalDate> preHmfArrivalDates) {
        this.dataPerSample = dataPerSample;
        this.preHmfArrivalDates = preHmfArrivalDates;
    }

    @Nullable
    public LocalDate arrivalDateForSample(@NotNull final String sample) {
        final LimsJsonData sampleData = dataPerSample.get(sample);
        final LocalDate arrivalDate;
        if (sampleData != null) {
            arrivalDate = getNullableDate(sampleData.arrivalDateString());
            if (arrivalDate == null) {
                LOGGER.warn("LIMS arrival date: " + sampleData.arrivalDateString() + " is not a valid date.");
            }
        } else {
            arrivalDate = preHmfArrivalDates.get(sample);
        }

        if (arrivalDate == null) {
            LOGGER.warn("Could not find arrival date for sample: " + sample + " in LIMS");
        }
        return arrivalDate;
    }

    @Nullable
    public LocalDate samplingDateForSample(@NotNull final String sample) {
        final LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            final LocalDate samplingDate = getNullableDate(sampleData.samplingDateString());
            if (samplingDate == null) {
                LOGGER.warn("LIMS sampling date: " + sampleData.samplingDateString() + " is not a valid date.");
            }
            return samplingDate;
        }
        LOGGER.warn("Could not find sampling date for sample: " + sample + " in LIMS");
        return null;
    }

    @Nullable
    public Double tumorPercentageForSample(@NotNull final String sample) {
        final LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            try {
                return Double.parseDouble(sampleData.tumorPercentage()) / 100D;
            } catch (final NumberFormatException e) {
                return null;
            }
        }
        LOGGER.warn("Could not find tumor percentage for sample: " + sample + " in LIMS");
        return null;
    }

    @NotNull
    public String labProceduresForSample(@NotNull final String sample) {
        final LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            return sampleData.labProcedures();
        }
        LOGGER.warn("Could not find lab SOP versions for sample: " + sample + " in LIMS");
        return "N/A";
    }

    @Nullable
    private static LocalDate getNullableDate(@NotNull final String dateString) {
        try {
            return LocalDate.parse(dateString, DATE_FORMATTER);
        } catch (DateTimeParseException e) {
            return null;
        }
    }
}
