package com.hartwig.hmftools.common.lims;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Lims {

    private static final Logger LOGGER = LogManager.getLogger(Lims.class);

    private Lims() {
    }

    @NotNull
    public static LimsModel buildModelFromCsv(@NotNull final String pathToCsv,
            @NotNull final DateTimeFormatter dateFormatter) throws IOException, EmptyFileException {
        final Map<String, LimsData> limsDataPerSample = Maps.newHashMap();
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines) {
            final String[] parts = line.split(",");
            final String sample = parts[0].trim();
            readLimsLine(sample, line, dateFormatter).ifPresent(
                    limsBiopsyData -> limsDataPerSample.put(sample, limsBiopsyData));
        }
        return new LimsModel(limsDataPerSample);
    }

    @NotNull
    private static Optional<LimsData> readLimsLine(@NotNull final String sample, @NotNull final String line,
            @NotNull final DateTimeFormatter dateFormatter) {
        final String[] parts = line.split(",");
        try {
            final String samplingDateField = !parts[1].trim().equals(Strings.EMPTY) ? parts[1].trim() : null;
            final String arrivalDateField = parts[2].trim();
            final LocalDate arrivalDate = LocalDate.parse(arrivalDateField, dateFormatter);
            final LocalDate samplingDate = getNullableDate(samplingDateField, dateFormatter);
            if (isReference(sample)) {
                if (parts.length == 4) {
                    LOGGER.warn("Sample " + sample + " has tumor percentage even though it is a reference sample.");
                }
                return Optional.of(new LimsBloodData(samplingDate, arrivalDate));
            }
            if (isTumor(sample)) {
                final Double tumorPercentage = parts.length == 4 ? readTumorPercentage(sample, parts[3]) : null;
                return Optional.of(new LimsTumorData(samplingDate, arrivalDate, tumorPercentage));
            }
            LOGGER.warn("Invalid sample name: " + sample + " in row: " + line);
            return Optional.empty();
        } catch (Exception e) {
            LOGGER.warn("Invalid row found in lims csv: " + line);
            return Optional.empty();
        }

    }

    @Nullable
    private static Double readTumorPercentage(@NotNull final String sample,
            @NotNull final String tumorPercentageField) {
        if (!tumorPercentageField.replaceAll("\"", Strings.EMPTY).trim().equals(Strings.EMPTY)) {
            try {
                return Double.valueOf(tumorPercentageField) / 100D;
            } catch (final NumberFormatException e) {
                LOGGER.warn(
                        "Could not parse tumor percentage from " + tumorPercentageField + " for sample: " + sample);
                return null;
            }
        }
        return null;
    }

    private static boolean isTumor(@NotNull final String sampleName) {
        return sampleName.length() >= 13 && sampleName.toLowerCase().charAt(12) == 't';
    }

    private static boolean isReference(@NotNull final String sampleName) {
        return sampleName.length() >= 13 && sampleName.toLowerCase().charAt(12) == 'r';
    }

    @Nullable
    private static LocalDate getNullableDate(@Nullable final String dateField,
            @NotNull final DateTimeFormatter dateFormatter) {
        if (dateField == null) {
            return null;
        }

        try {
            return LocalDate.parse(dateField, dateFormatter);
        } catch (DateTimeParseException e) {
            return null;
        }
    }

    @NotNull
    public static LimsModel buildEmptyModel() {
        return new LimsModel(Maps.newHashMap());
    }
}
