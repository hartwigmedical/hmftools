package com.hartwig.hmftools.common.lims;

import java.io.File;
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
import com.google.common.collect.Sets;
import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonSyntaxException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LimsFactory {

    private static final Logger LOGGER = LogManager.getLogger(LimsFactory.class);

    private static final String LIMS_JSON_FILE = "lims.json";
    private static final String PRE_LIMS_ARRIVAL_DATES_FILE = "pre_lims_arrival_dates.csv";
    private static final String SAMPLES_WITHOUT_SAMPLING_DATE_FILE = "samples_without_sampling_date.csv";

    private LimsFactory() {
    }

    @NotNull
    public static Lims fromLimsDirectory(@NotNull final String limsDirectory) throws IOException {
        final String limsJsonPath = limsDirectory + File.separator + LIMS_JSON_FILE;
        Map<String, LimsJsonSampleData> dataPerSample = readLimsJsonSamples(limsJsonPath);
        Map<String, LimsJsonSubmissionData> dataPerSubmission = readLimsJsonSubmissions(limsJsonPath);

        Map<String, LocalDate> preLIMSArrivalDates =
                readPreLIMSArrivalDateCsv(limsDirectory + File.separator + PRE_LIMS_ARRIVAL_DATES_FILE);
        Set<String> samplesWithoutSamplingDate =
                readSamplesWithoutSamplingDateCsv(limsDirectory + File.separator + SAMPLES_WITHOUT_SAMPLING_DATE_FILE);
        return new Lims(dataPerSample, dataPerSubmission, preLIMSArrivalDates, samplesWithoutSamplingDate);
    }

    @NotNull
    public static Lims empty() {
        return new Lims(Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap(), Sets.newHashSet());
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsJsonSubmissionData> readLimsJsonSubmissions(@NotNull final String limsJsonPath)
            throws FileNotFoundException {
        final Gson gson = LimsGsonAdapter.buildSubmissionGson();
        final JsonObject jsonObject = new JsonParser().parse(new FileReader(limsJsonPath)).getAsJsonObject();
        final Set<Map.Entry<String, JsonElement>> jsonSubmissions = jsonObject.getAsJsonObject("submissions").entrySet();
        final Map<String, LimsJsonSubmissionData> limsDataPerSubmission = Maps.newHashMap();

        jsonSubmissions.forEach(jsonSubmission -> {
            final JsonObject jsonSampleObject = jsonSubmission.getValue().getAsJsonObject();
            final String projectType = jsonSampleObject.get("project_type").getAsString();
            if (projectType.contains("CORE")) {
                try {
                    final LimsJsonSubmissionData limsJsonSubmissionData =
                            gson.fromJson(jsonSubmission.getValue(), LimsJsonSubmissionData.class);
                    limsDataPerSubmission.put(limsJsonSubmissionData.submission(), limsJsonSubmissionData);
                } catch (JsonSyntaxException e) {
                    LOGGER.warn("Could not convert json element to LimsJsonSubmissionData: " + jsonSubmission.getValue() + " - message:" + e
                            .getMessage());
                }
            }
        });

        return limsDataPerSubmission;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsJsonSampleData> readLimsJsonSamples(@NotNull final String limsJsonPath) throws FileNotFoundException {
        final Gson gson = LimsGsonAdapter.buildSampleGson();
        final JsonObject jsonObject = new JsonParser().parse(new FileReader(limsJsonPath)).getAsJsonObject();
        final Set<Map.Entry<String, JsonElement>> jsonSamples = jsonObject.getAsJsonObject("samples").entrySet();
        final Map<String, LimsJsonSampleData> limsDataPerSample = Maps.newHashMap();

        jsonSamples.forEach(jsonSample -> {
            final JsonObject jsonSampleObject = jsonSample.getValue().getAsJsonObject();
            final String analysisType = jsonSampleObject.get("analysis_type").getAsString();
            final String label = jsonSampleObject.get("label").getAsString();

            // Filter on somatic to get rid of RNA samples, see also DEV-252
            // We are not interested in research labeled samples.
            if (analysisType != null && analysisType.toLowerCase().contains("somatic") && !label.equalsIgnoreCase("research")) {
                try {
                    final LimsJsonSampleData limsJsonSampleData = gson.fromJson(jsonSample.getValue(), LimsJsonSampleData.class);
                    limsDataPerSample.put(limsJsonSampleData.sampleId(), limsJsonSampleData);
                } catch (JsonSyntaxException e) {
                    LOGGER.warn("Could not convert json element to LimsJsonSampleData: " + jsonSample.getValue() + " - message:"
                            + e.getMessage());
                }
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
                    arrivalDate = LocalDate.parse(arrivalDateString, LimsConstants.DATE_FORMATTER);
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

    @NotNull
    @VisibleForTesting
    static Set<String> readSamplesWithoutSamplingDateCsv(@NotNull String samplesWithoutSamplingDate) throws IOException {
        return Files.lines(Paths.get(samplesWithoutSamplingDate)).filter(s -> !s.isEmpty()).collect(Collectors.toSet());
    }
}
