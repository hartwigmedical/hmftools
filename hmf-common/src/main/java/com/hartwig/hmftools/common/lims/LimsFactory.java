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
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModelFactory;
import com.hartwig.hmftools.common.lims.hospital.HospitalModel;
import com.hartwig.hmftools.common.lims.hospital.HospitalModelFactory;
import com.hartwig.hmftools.common.lims.hospital.ImmutableHospitalModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LimsFactory {

    private static final Logger LOGGER = LogManager.getLogger(LimsFactory.class);

    private static final String LIMS_JSON_FILE = "lims.json";

    private static final String PRE_LIMS_ARRIVAL_DATES_TSV = "pre_lims_arrival_dates.tsv";
    private static final String SAMPLES_WITHOUT_SAMPLING_DATE_TSV = "samples_without_sampling_date.tsv";
    private static final String LIMS_SHALLOW_SEQ_TSV = "shallow_seq_purity.tsv";
    private static final String PATIENT_BLACKLIST_TSV = "patient_blacklist.tsv";
    private static final String PATIENTS_WITHOUT_CURATED_PRIMARY_TUMOR = "patients_without_curated_primary_tumor.tsv";
    private static final String COHORT_CONFIG_TSV = "cohort_config.tsv";

    private static final String FIELD_SEPARATOR = "\t";

    private LimsFactory() {
    }

    @NotNull
    public static Lims fromLimsDirectory(@NotNull String limsDirectory) throws IOException {
        String limsJsonPath = limsDirectory + File.separator + LIMS_JSON_FILE;
        Map<String, LimsJsonSampleData> dataPerSampleBarcode = readLimsJsonSamples(limsJsonPath);
        Map<String, LimsJsonSubmissionData> dataPerSubmission = readLimsJsonSubmissions(limsJsonPath);
        Map<String, LimsShallowSeqData> shallowSeqPerSampleBarcode =
                readLimsShallowSeqTsv(limsDirectory + File.separator + LIMS_SHALLOW_SEQ_TSV);

        Map<String, LocalDate> preLimsArrivalDates = readPreLimsArrivalDateTsv(limsDirectory + File.separator + PRE_LIMS_ARRIVAL_DATES_TSV);
        Set<String> samplesWithoutSamplingDate = readSingleColumnTsv(limsDirectory + File.separator + SAMPLES_WITHOUT_SAMPLING_DATE_TSV);
        Set<String> blacklistedPatients = readSingleColumnTsv(limsDirectory + File.separator + PATIENT_BLACKLIST_TSV);
        Set<String> patientsWithoutCuratedPrimaryTumor =
                readSingleColumnTsv(limsDirectory + File.separator + PATIENTS_WITHOUT_CURATED_PRIMARY_TUMOR);

        HospitalModel hospitalModel = HospitalModelFactory.fromLimsDirectory(limsDirectory);

        LimsCohortModel cohortModel = LimsCohortModelFactory.read(limsDirectory + File.separator + COHORT_CONFIG_TSV);

        return new Lims(cohortModel,
                hospitalModel,
                dataPerSampleBarcode,
                dataPerSubmission,
                shallowSeqPerSampleBarcode,
                preLimsArrivalDates,
                samplesWithoutSamplingDate,
                patientsWithoutCuratedPrimaryTumor,
                blacklistedPatients);
    }

    @NotNull
    public static Lims empty() {
        LimsCohortModel alwaysDisabledCohortModel = new LimsCohortModel() {
            @NotNull
            @Override
            protected Map<String, LimsCohortConfig> limsCohortMap() {
                return Maps.newHashMap();
            }

            @Nullable
            @Override
            public LimsCohortConfig queryCohortData(@Nullable final String cohortString, @NotNull final String sampleId) {
                return ImmutableLimsCohortConfig.builder()
                        .cohortId(sampleId)
                        .sampleContainsHospitalCenterId(false)
                        .reportGermline(false)
                        .reportGermlineFlag(false)
                        .reportConclusion(false)
                        .reportViral(false)
                        .reportPeach(false)
                        .requireHospitalId(false)
                        .requireHospitalPAId(false)
                        .requireHospitalPersonsStudy(false)
                        .requireHospitalPersonsRequester(false)
                        .requireAdditionalInformationForSidePanel(false)
                        .build();
            }
        };

        return new Lims(alwaysDisabledCohortModel,
                ImmutableHospitalModel.builder().build(),
                Maps.newHashMap(),
                Maps.newHashMap(),
                Maps.newHashMap(),
                Maps.newHashMap(),
                Sets.newHashSet(),
                Sets.newHashSet(),
                Sets.newHashSet());
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsJsonSampleData> readLimsJsonSamples(@NotNull String limsJsonPath) throws FileNotFoundException {
        Gson gson = LimsGsonAdapter.buildSampleGson();
        JsonObject jsonObject = new JsonParser().parse(new FileReader(limsJsonPath)).getAsJsonObject();
        Set<Map.Entry<String, JsonElement>> jsonSamples = jsonObject.getAsJsonObject("samples").entrySet();

        Map<String, LimsJsonSampleData> limsDataPerSampleBarcode = Maps.newHashMap();

        jsonSamples.forEach(jsonSample -> {
            String barcode = jsonSample.getKey();
            JsonObject jsonSampleObject = jsonSample.getValue().getAsJsonObject();
            String analysisType = jsonSampleObject.get("analysis_type").getAsString();
            String label = jsonSampleObject.get("label").getAsString();

            // DEV-252 - Filter on somatic to get rid of RNA samples
            // Also, we are not interested in research-labeled samples.
            if (analysisType != null && analysisType.toLowerCase().contains("somatic") && !label.equalsIgnoreCase("research")) {
                try {
                    LimsJsonSampleData limsJsonSampleData = gson.fromJson(jsonSample.getValue(), LimsJsonSampleData.class);
                    if (limsDataPerSampleBarcode.containsKey(barcode)) {
                        LOGGER.warn("LIMS contains duplicate entries for {} ({})", barcode, limsJsonSampleData.sampleId());
                    }
                    limsDataPerSampleBarcode.put(barcode, limsJsonSampleData);
                } catch (JsonSyntaxException exception) {
                    LOGGER.warn("Could not convert JSON element to LimsJsonSampleData: {} - message: {}",
                            jsonSample.getValue(),
                            exception.getMessage());
                }
            }
        });
        //Make COLO patient ID available for use cases
        limsDataPerSampleBarcode.put("COLO829V003TVAL", createLimsSampleDataForCOLO());

        return limsDataPerSampleBarcode;
    }

    @NotNull
    private static LimsJsonSampleData createLimsSampleDataForCOLO() {
        return ImmutableLimsJsonSampleData.builder()
                .sampleId(Strings.EMPTY)
                .patientId("COLO829")
                .tumorBarcode("COLO829V003TVAL")
                .refBarcode(Strings.EMPTY)
                .arrivalDate(Strings.EMPTY)
                .dnaConcentration(Strings.EMPTY)
                .primaryTumor(Strings.EMPTY)
                .labSopVersions(Strings.EMPTY)
                .submission(Strings.EMPTY)
                .germlineReportingLevel(Strings.EMPTY)
                .reportGermlineVariants(false)
                .shallowSeq(false)
                .reportViralInsertions(false)
                .cohort(Strings.EMPTY)
                .analysisType(Strings.EMPTY).build();
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsJsonSubmissionData> readLimsJsonSubmissions(@NotNull String limsJsonPath) throws FileNotFoundException {
        Gson gson = LimsGsonAdapter.buildSubmissionGson();
        JsonObject jsonObject = new JsonParser().parse(new FileReader(limsJsonPath)).getAsJsonObject();
        Set<Map.Entry<String, JsonElement>> jsonSubmissions = jsonObject.getAsJsonObject("submissions").entrySet();

        Map<String, LimsJsonSubmissionData> limsDataPerSubmission = Maps.newHashMap();

        jsonSubmissions.forEach(jsonSubmission -> {
            JsonObject jsonSampleObject = jsonSubmission.getValue().getAsJsonObject();
            String projectType = jsonSampleObject.get("project_type").getAsString();
            // We only need submission data for CORE projects
            if (projectType.contains("CORE")) {
                try {
                    LimsJsonSubmissionData limsJsonSubmissionData = gson.fromJson(jsonSubmission.getValue(), LimsJsonSubmissionData.class);
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
    static Map<String, LocalDate> readPreLimsArrivalDateTsv(@NotNull String preLimsArrivalDatesTsv) throws IOException {
        Map<String, LocalDate> arrivalDatesPerSampleId = Maps.newHashMap();
        List<String> lines = Files.lines(Paths.get(preLimsArrivalDatesTsv)).collect(Collectors.toList());
        for (String line : lines) {
            String[] parts = line.split(FIELD_SEPARATOR);

            if (parts.length == 2) {
                String sampleId = parts[0].trim();
                String arrivalDateString = parts[1].trim();
                LocalDate arrivalDate;
                try {
                    arrivalDate = LocalDate.parse(arrivalDateString, LimsConstants.DATE_FORMATTER);
                } catch (DateTimeParseException exc) {
                    LOGGER.warn("Could not parse date in pre-HMF arrival date csv: {}", arrivalDateString);
                    arrivalDate = null;
                }
                arrivalDatesPerSampleId.put(sampleId, arrivalDate);
            } else {
                LOGGER.warn("Invalid line in pre-HMF arrival date csv: {}", line);
            }
        }
        return arrivalDatesPerSampleId;
    }

    @NotNull
    @VisibleForTesting
    static Map<String, LimsShallowSeqData> readLimsShallowSeqTsv(@NotNull String shallowSeqTsv) throws IOException {
        Map<String, LimsShallowSeqData> shallowSeqPerSampleBarcode = Maps.newHashMap();
        List<String> lines = Files.lines(Paths.get(shallowSeqTsv)).collect(Collectors.toList());

        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, 5);
            String sampleBarcode = parts[0];
            if (parts.length == 5) {
                shallowSeqPerSampleBarcode.put(sampleBarcode,
                        ImmutableLimsShallowSeqData.builder()
                                .sampleBarcode(sampleBarcode)
                                .sampleId(parts[1])
                                .purityShallowSeq(parts[2])
                                .hasReliableQuality(Boolean.parseBoolean(parts[3]))
                                .hasReliablePurity(Boolean.parseBoolean(parts[4]))
                                .build());
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in shallow seq csv: {}", line);
            }
        }
        return shallowSeqPerSampleBarcode;
    }

    @NotNull
    @VisibleForTesting
    static Set<String> readSingleColumnTsv(@NotNull String tsv) throws IOException {
        return Files.lines(Paths.get(tsv)).filter(s -> !s.isEmpty()).collect(Collectors.toSet());
    }
}
