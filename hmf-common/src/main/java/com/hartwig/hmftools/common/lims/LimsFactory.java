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
    private static final String PRE_LIMS_ARRIVAL_DATES_TSV = "pre_lims_arrival_dates_pilot.tsv";
    private static final String SAMPLES_WITHOUT_SAMPLING_DATE_TSV = "samples_without_sampling_date_pilot.tsv";
    private static final String LIMS_SHALLOW_SEQ_TSV = "shallow_seq_purity_pilot.tsv";

    private static final String FIELD_SEPARATOR = "\t";

    private LimsFactory() {
    }

    @NotNull
    public static Lims fromLimsDirectory(@NotNull final String limsDirectory) throws IOException {
        String limsJsonPath = limsDirectory + File.separator + LIMS_JSON_FILE;
        Map<String, LimsJsonSampleData> dataPerSampleBarcode = readLimsJsonSamples(limsJsonPath);
        Map<String, LimsJsonSubmissionData> dataPerSubmission = readLimsJsonSubmissions(limsJsonPath);

        Map<String, LocalDate> preLIMSArrivalDates = readPreLIMSArrivalDateTsv(limsDirectory + File.separator + PRE_LIMS_ARRIVAL_DATES_TSV);
        Set<String> samplesWithoutSamplingDate =
                readSamplesWithoutSamplingDateTsv(limsDirectory + File.separator + SAMPLES_WITHOUT_SAMPLING_DATE_TSV);
        Map<String, LimsShallowSeqData> shallowSeqPerSampleBarcode =
                readLimsShallowSeq(limsDirectory + File.separator + LIMS_SHALLOW_SEQ_TSV);

        return new Lims(dataPerSampleBarcode,
                dataPerSubmission,
                preLIMSArrivalDates,
                samplesWithoutSamplingDate,
                shallowSeqPerSampleBarcode);
    }

    @NotNull
    public static Lims empty() {
        return new Lims(Maps.newHashMap(), Maps.newHashMap(), Maps.newHashMap(), Sets.newHashSet(), Maps.newHashMap());
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsShallowSeqData> readLimsShallowSeq(@NotNull final String shallowSeqTsv) throws IOException {
        final Map<String, LimsShallowSeqData> shallowSeqPerSampleBarcode = Maps.newHashMap();
        final List<String> lines =
                com.hartwig.hmftools.common.utils.io.reader.FileReader.build().readLines(new File(shallowSeqTsv).toPath());
        for (final String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, 5);
            if (parts.length == 5) {
                shallowSeqPerSampleBarcode.put(parts[0],
                        ImmutableLimsShallowSeqData.of(parts[0],
                                parts[1],
                                parts[2],
                                Boolean.parseBoolean(parts[3]),
                                Boolean.parseBoolean(parts[4])));
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in shallow seq csv: {}", line);
            }
        }
        return shallowSeqPerSampleBarcode;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsJsonSubmissionData> readLimsJsonSubmissions(@NotNull final String limsJsonPath) throws FileNotFoundException {
        final Gson gson = LimsGsonAdapter.buildSubmissionGson();
        final JsonObject jsonObject = new JsonParser().parse(new FileReader(limsJsonPath)).getAsJsonObject();
        final Set<Map.Entry<String, JsonElement>> jsonSubmissions = jsonObject.getAsJsonObject("submissions").entrySet();
        final Map<String, LimsJsonSubmissionData> limsDataPerSubmission = Maps.newHashMap();

        jsonSubmissions.forEach(jsonSubmission -> {
            final JsonObject jsonSampleObject = jsonSubmission.getValue().getAsJsonObject();
            final String projectType = jsonSampleObject.get("project_type").getAsString();
            // We only need submission data for CORE projects for now
            if (projectType.contains("CORE")) {
                try {
                    final LimsJsonSubmissionData limsJsonSubmissionData =
                            gson.fromJson(jsonSubmission.getValue(), LimsJsonSubmissionData.class);
                    limsDataPerSubmission.put(limsJsonSubmissionData.submission(), limsJsonSubmissionData);
                } catch (JsonSyntaxException exception) {
                    LOGGER.warn("Could not convert JSON element to LimsJsonSubmissionData: {} - message: {}",
                            jsonSubmission.getValue(),
                            exception.getMessage());
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
            final String barcode = jsonSample.getKey();
            final JsonObject jsonSampleObject = jsonSample.getValue().getAsJsonObject();
            final String analysisType = jsonSampleObject.get("analysis_type").getAsString();
            final String label = jsonSampleObject.get("label").getAsString();

            // DEV-252 - Filter on somatic to get rid of RNA samples
            // Also, we are not interested in research-labeled samples.
            if (analysisType != null && analysisType.toLowerCase().contains("somatic") && !label.equalsIgnoreCase("research")) {
                try {
                    final LimsJsonSampleData limsJsonSampleData = gson.fromJson(jsonSample.getValue(), LimsJsonSampleData.class);
                    if (limsDataPerSample.containsKey(barcode)) {
                        LOGGER.warn("LIMS contains duplicate entries for {} ({})", barcode, limsJsonSampleData.sampleId());
                    }
                    limsDataPerSample.put(barcode, limsJsonSampleData);
                } catch (JsonSyntaxException exception) {
                    LOGGER.warn("Could not convert JSON element to LimsJsonSampleData: {} - message: {}",
                            jsonSample.getValue(),
                            exception.getMessage());
                }
            }
        });

        return limsDataPerSample;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LocalDate> readPreLIMSArrivalDateTsv(@NotNull final String preLIMSArrivalDatesTsvPath) throws IOException {
        final Map<String, LocalDate> arrivalDatesPerSample = Maps.newHashMap();
        final List<String> lines = Files.lines(Paths.get(preLIMSArrivalDatesTsvPath)).collect(Collectors.toList());
        for (final String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR);

            if (parts.length == 2) {
                final String sample = parts[0].trim();
                final String arrivalDateString = parts[1].trim();
                LocalDate arrivalDate;
                try {
                    arrivalDate = LocalDate.parse(arrivalDateString, LimsConstants.DATE_FORMATTER);
                } catch (DateTimeParseException exc) {
                    LOGGER.warn("Could not parse date in pre-HMF arrival date csv: {}", arrivalDateString);
                    arrivalDate = null;
                }
                arrivalDatesPerSample.put(sample, arrivalDate);
            } else {
                LOGGER.warn("Invalid line in pre-HMF arrival date csv: {}", line);
            }
        }
        return arrivalDatesPerSample;
    }

    @NotNull
    @VisibleForTesting
    static Set<String> readSamplesWithoutSamplingDateTsv(@NotNull String samplesWithoutSamplingDate) throws IOException {
        return Files.lines(Paths.get(samplesWithoutSamplingDate)).filter(s -> !s.isEmpty()).collect(Collectors.toSet());
    }
}
