package com.hartwig.hmftools.common.lims;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.LocalDate;
import java.time.format.DateTimeParseException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LimsFactory {

    private static final Logger LOGGER = LogManager.getLogger(Lims.class);

    private LimsFactory() {
    }

    @NotNull
    public static Lims fromLimsJson(@NotNull final String limsJsonPath) throws FileNotFoundException {
        return new Lims(readLimsJson(limsJsonPath), Maps.newHashMap());
    }

    @NotNull
    public static Lims fromLimsJsonWithPreLIMSArrivalDates(@NotNull final String limsJsonPath,
            @NotNull final String preLIMSArrivalDatesCsvPath) throws IOException {
        return new Lims(readLimsJson(limsJsonPath), readPreLIMSArrivalDateCsv(preLIMSArrivalDatesCsvPath));
    }

    @NotNull
    public static Lims empty() {
        return new Lims(Maps.newHashMap(), Maps.newHashMap());
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsJsonData> readLimsJson(@NotNull final String limsJsonPath) throws FileNotFoundException {
        final Gson gson = LimsGsonAdapter.buildGson();
        final JsonObject jsonObject = new JsonParser().parse(new FileReader(limsJsonPath)).getAsJsonObject();
        final Set<Map.Entry<String, JsonElement>> jsonSamples = jsonObject.getAsJsonObject("samples").entrySet();
        final Map<String, LimsJsonData> limsDataPerSample = Maps.newHashMap();

        jsonSamples.forEach(jsonSample -> {
            final JsonObject jsonSampleObject = jsonSample.getValue().getAsJsonObject();
            final String sampleLabel = jsonSampleObject.get("label").getAsString();
            final String analysisType = jsonSampleObject.get("analysis_type").getAsString();

            if ((sampleLabel.equals("CPCT") || sampleLabel.equals("DRUP")) && analysisType != null && analysisType.toLowerCase()
                    .contains("somatic")) {
                final LimsJsonData limsJsonData = gson.fromJson(jsonSample.getValue(), LimsJsonData.class);
                limsDataPerSample.put(limsJsonData.sampleId(), limsJsonData);
            }
        });

        return limsDataPerSample;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LocalDate> readPreLIMSArrivalDateCsv(@NotNull final String preLIMSArrivalDatesCsvPath) throws IOException {
        final Map<String, LocalDate> arrivalDatesPerSample = Maps.newHashMap();
        final List<String> lines = Files.lines(Paths.get(preLIMSArrivalDatesCsvPath)).collect(Collectors.toList());
        for (final String line : lines) {
            final String[] parts = line.split(",");

            if (parts.length == 2) {
                final String sample = parts[0].trim();
                final String arrivalDateString = parts[1].trim();
                LocalDate arrivalDate;
                try {
                    arrivalDate = LocalDate.parse(arrivalDateString, Lims.DATE_FORMATTER);
                } catch (DateTimeParseException exc) {
                    LOGGER.warn("Could not parse date in pre-HMF arrival date csv: " + arrivalDateString);
                    arrivalDate = null;
                }
                arrivalDatesPerSample.put(sample, arrivalDate);
            } else {
                LOGGER.warn("Invalid line in pre-HMF arrival date csv: " + line);
            }
        }
        return arrivalDatesPerSample;
    }
}
