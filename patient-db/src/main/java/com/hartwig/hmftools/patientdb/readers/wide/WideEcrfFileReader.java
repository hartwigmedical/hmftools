package com.hartwig.hmftools.patientdb.readers.wide;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class WideEcrfFileReader {

    private static final Logger LOGGER = LogManager.getLogger(WideEcrfFileReader.class);
    private static final String FIELD_SEPARATOR = ";";

    private static final int PRE_TREATMENT_DATA_COUNT = 7;
    private static final int PRE_TREATMENT_DATA_PATIENT_ID = 0;
    private static final int PRE_TREATMENT_DATA_PREVIOUS_THERAPY = 1;
    private static final int PRE_TREATMENT_DATA_DRUG1 = 2;
    private static final int PRE_TREATMENT_DATA_DRUG2 = 3;
    private static final int PRE_TREATMENT_DATA_DRUG3 = 4;
    private static final int PRE_TREATMENT_DATA_DRUG4 = 5;
    private static final int PRE_TREATMENT_DATA_DATE_LAST_SYSTEMIC_THERAPY = 6;

    private static final int BIOPSY_DATA_COUNT = 6;
    private static final int BIOPSY_DATA_PATIENT_ID = 0;
    private static final int BIOPSY_DATA_WIDE_ID = 1;
    private static final int BIOPSY_DATA_DATA_AVAILABLE = 2;
    private static final int BIOPSY_DATA_TISSUE_ID = 3;
    private static final int BIOPSY_DATA_DATE = 4;
    private static final int BIOPSY_DATA_WGS_SUCCESSFUL = 5;

    private static final int TREATMENT_DATA_COUNT = 5;
    private static final int TREATMENT_DATA_ID = 0;
    private static final int TREATMENT_DATA_CODE = 1;
    private static final int TREATMENT_DATA_DRUG = 2;
    private static final int TREATMENT_DATA_START_DATE = 3;
    private static final int TREATMENT_DATA_END_DATE = 4;

    private static final int RESPONSE_DATA_COUNT = 8;
    private static final int RESPONSE_DATA_PATIENT_ID = 0;
    private static final int RESPONSE_DATA_TIME_POINT = 1;
    private static final int RESPONSE_DATA_DATE = 2;
    private static final int RESPONSE_DATA_RECIST_NOT_DONE = 3;
    private static final int RESPONSE_DATA_RESPONSE_ACCORDING_RECIST = 4;
    private static final int RESPONSE_DATA_CLINICAL_DECISION = 5;
    private static final int RESPONSE_DATA_REASON_STOP_TREATMENT = 6;
    private static final int RESPONSE_DATA_REASON_STOP_TREATMENT_OTHER = 7;

    private WideEcrfFileReader() {
    }

    @NotNull
    public static List<WideTreatmentData> readTreatmentData(@NotNull String pathToCsv) throws IOException {
        List<WideTreatmentData> wideTreatments = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, TREATMENT_DATA_COUNT);
            if (parts.length == TREATMENT_DATA_COUNT) {
                WideTreatmentData wideTreatment = ImmutableWideTreatmentData.of(parts[TREATMENT_DATA_ID].replace("-", ""),
                        parts[TREATMENT_DATA_CODE],
                        parts[TREATMENT_DATA_DRUG],
                        parts[TREATMENT_DATA_START_DATE],
                        parts[TREATMENT_DATA_END_DATE]);

                wideTreatments.add(wideTreatment);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in WIDE treatment csv: {}", line);
            }
        }
        LOGGER.info(wideTreatments);
        return wideTreatments;
    }

    @NotNull
    public static List<WidePreTreatmentData> readPreTreatmentData(@NotNull String pathToCsv) throws IOException {
        List<WidePreTreatmentData> widePreTreatments = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, PRE_TREATMENT_DATA_COUNT);
            if (parts.length == PRE_TREATMENT_DATA_COUNT) {
                WidePreTreatmentData widePreTreatmentData =
                        ImmutableWidePreTreatmentData.of(parts[PRE_TREATMENT_DATA_PATIENT_ID].replace("-", ""),
                                parts[PRE_TREATMENT_DATA_PREVIOUS_THERAPY],
                                parts[PRE_TREATMENT_DATA_DRUG1],
                                parts[PRE_TREATMENT_DATA_DRUG2],
                                parts[PRE_TREATMENT_DATA_DRUG3],
                                parts[PRE_TREATMENT_DATA_DRUG4],
                                parts[PRE_TREATMENT_DATA_DATE_LAST_SYSTEMIC_THERAPY]);

                widePreTreatments.add(widePreTreatmentData);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in WIDE pre treatment data csv: {}", line);
            }
        }
        LOGGER.info(widePreTreatments);
        return widePreTreatments;
    }

    @NotNull
    public static List<WideBiopsyData> readBiopsyData(@NotNull String pathToCsv) throws IOException {
        List<WideBiopsyData> wideBiopsies = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, BIOPSY_DATA_COUNT);
            if (parts.length == BIOPSY_DATA_COUNT) {
                WideBiopsyData wideBiopsyData = ImmutableWideBiopsyData.of(parts[BIOPSY_DATA_PATIENT_ID],
                        parts[BIOPSY_DATA_WIDE_ID].replace("-", ""),
                        parts[BIOPSY_DATA_DATA_AVAILABLE],
                        parts[BIOPSY_DATA_TISSUE_ID],
                        parts[BIOPSY_DATA_DATE],
                        parts[BIOPSY_DATA_WGS_SUCCESSFUL]);

                wideBiopsies.add(wideBiopsyData);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in WIDE biopsy date csv: {}", line);
            }
        }
        LOGGER.info(wideBiopsies);
        return wideBiopsies;
    }

    @NotNull
    public static List<WideResponseData> readResponseData(@NotNull String pathToCsv) throws IOException {
        List<WideResponseData> wideResponses = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_SEPARATOR, RESPONSE_DATA_COUNT);
            if (parts.length == RESPONSE_DATA_COUNT) {
                WideResponseData wideResponseData = ImmutableWideResponseData.of(parts[RESPONSE_DATA_PATIENT_ID].replace("-", ""),
                        parts[RESPONSE_DATA_TIME_POINT],
                        parts[RESPONSE_DATA_DATE],
                        parts[RESPONSE_DATA_RECIST_NOT_DONE],
                        parts[RESPONSE_DATA_RESPONSE_ACCORDING_RECIST],
                        parts[RESPONSE_DATA_CLINICAL_DECISION],
                        parts[RESPONSE_DATA_REASON_STOP_TREATMENT],
                        parts[RESPONSE_DATA_REASON_STOP_TREATMENT_OTHER]);

                wideResponses.add(wideResponseData);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in WIDE biopsy date csv: {}", line);
            }
        }
        LOGGER.info(wideResponses);
        return wideResponses;
    }
}
