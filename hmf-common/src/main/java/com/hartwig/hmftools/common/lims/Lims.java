package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;
import java.time.format.DateTimeParseException;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Lims {

    private static final Logger LOGGER = LogManager.getLogger(Lims.class);

    @NotNull
    private final Map<String, LimsJsonData> dataPerSample;
    @NotNull
    private final Map<String, LocalDate> preLimsArrivalDates;

    Lims(@NotNull final Map<String, LimsJsonData> dataPerSample, @NotNull final Map<String, LocalDate> preLimsArrivalDates) {
        this.dataPerSample = dataPerSample;
        this.preLimsArrivalDates = preLimsArrivalDates;
    }

    public int sampleCount() {
        return dataPerSample.size();
    }

    @Nullable
    public LocalDate arrivalDateForSample(@NotNull final String sample) {
        LimsJsonData sampleData = dataPerSample.get(sample);
        final LocalDate arrivalDate;
        if (sampleData != null) {
            arrivalDate = getNullableDate(sampleData.arrivalDateString());
            if (arrivalDate == null) {
                LOGGER.warn("LIMS arrival date for " + sample + ": " + sampleData.arrivalDateString() + " is not a valid date.");
            }
        } else {
            arrivalDate = preLimsArrivalDates.get(sample);
        }

        if (arrivalDate == null) {
            LOGGER.warn("Could not find arrival date for sample: " + sample + " in LIMS");
        }
        return arrivalDate;
    }

    @Nullable
    public LocalDate samplingDateForSample(@NotNull final String sample) {
        LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            final String samplingDateString = sampleData.samplingDateString();
            final LocalDate samplingDate = getNullableDate(samplingDateString);
            if (samplingDate == null && samplingDateString != null && !samplingDateString.equalsIgnoreCase("na")) {
                LOGGER.warn("LIMS sampling date for " + sample + ": " + sampleData.samplingDateString() + " is not a valid date.");
            }
            return samplingDate;
        }
        return null;
    }

    @Nullable
    public Integer dnaNanogramsForSample(@NotNull String sample) {
        LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            try {
                // KODU: LIMS stores the amount of nanograms per micro liter.
                return (int) Math.round(Double.parseDouble(sampleData.dnaConcentration()) * LimsConstants.DNA_MICRO_LITERS);
            } catch (final NumberFormatException e) {
                return null;
            }
        }
        return null;
    }

    @Nullable
    public Double tumorPercentageForSample(@NotNull final String sample) {
        LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            String tumorPercentageString = sampleData.tumorPercentageString();
            if (tumorPercentageString == null) {
                return null;
            }
            try {
                return Double.parseDouble(tumorPercentageString) / 100D;
            } catch (final NumberFormatException e) {
                return null;
            }
        }
        return null;
    }

    @NotNull
    public String primaryTumorForSample(@NotNull final String sample) {
        LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            return sampleData.primaryTumor();
        }
        // KODU: No warning raised since initially this information was not tracked so this will be missing for early samples.
        return "N/A";
    }

    @NotNull
    public String labProceduresForSample(@NotNull final String sample) {
        LimsJsonData sampleData = dataPerSample.get(sample);
        if (sampleData != null) {
            return sampleData.labProcedures();
        }
        LOGGER.warn("Could not find lab SOP versions for sample: " + sample + " in LIMS");
        return "N/A";
    }

    @Nullable
    private static LocalDate getNullableDate(@Nullable final String dateString) {
        if (dateString == null) {
            return null;
        }

        try {
            return LocalDate.parse(dateString, LimsConstants.DATE_FORMATTER);
        } catch (DateTimeParseException e) {
            return null;
        }
    }
}
