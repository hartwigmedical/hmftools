package com.hartwig.hmftools.patientdb.clinical.readers.wide;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class WideEcrfFileReader {

    private static final Logger LOGGER = LogManager.getLogger(WideEcrfFileReader.class);
    private static final char COMMA_DELIMITER = ',';
    private static final char SEMICOLON_DELIMITER = ';';

    private static final int FIVE_DAYS_DATA_COUNT = 15;
    private static final int FIVE_DAYS_DATA_PATIENT_ID = 1;
    private static final int FIVE_DAYS_DATA_BIRTH_YEAR = 2;
    private static final int FIVE_DAYS_DATA_GENDER = 3;
    private static final int FIVE_DAYS_DATA_INFORMED_CONSENT_DATE = 4;
    private static final int FIVE_DAYS_DATA_IS_AVAILABLE = 6;
    private static final int FIVE_DAYS_DATA_BIOPSY_DATE = 7;
    private static final int FIVE_DAYS_DATA_BIOPSY_SITE = 8;
    private static final int FIVE_DAYS_DATA_SAMPLE_TISSUE = 9;
    private static final int FIVE_DAYS_DATA_SAMPLE_TYPE = 10;
    private static final int FIVE_DAYS_DATA_STUDY_CODE = 11;
    private static final int FIVE_DAYS_DATA_OTHER_TRIALS = 12;
    private static final int FIVE_DAYS_DATA_OTHER_TRIAL_CODES = 13;
    private static final int FIVE_DAYS_DATA_OTHER_TRIAL_START_DATES = 14;

    private static final int PRE_AVL_TREATMENT_DATA_COUNT = 9;
    private static final int PRE_AVL_TREATMENT_DATA_PATIENT_ID = 0;
    private static final int PRE_AVL_TREATMENT_DATA_HAS_PREVIOUS_THERAPY = 2;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG1 = 3;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG2 = 4;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG3 = 5;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG4 = 6;
    private static final int PRE_AVL_TREATMENT_DATA_LAST_SYSTEMIC_THERAPY_DATE = 7;

    private static final int BIOPSY_DATA_COUNT = 6;
    private static final int BIOPSY_DATA_PATIENT_ID = 1;
    private static final int BIOPSY_DATA_PATHOLOGY_SAMPLE_ID = 3;
    private static final int BIOPSY_DATA_DATE = 4;
    private static final int BIOPSY_DATA_HAS_RECEIVED_SUCCESSFUL_REPORT = 5;

    private static final int AVL_TREATMENT_DATA_COUNT = 6;
    private static final int AVL_TREATMENT_DATA_PATIENT_ID = 0;
    private static final int AVL_TREATMENT_DATA_DRUG_CODE = 2;
    private static final int AVL_TREATMENT_DATA_DRUG = 3;
    private static final int AVL_TREATMENT_DATA_START_DATE = 4;
    private static final int AVL_TREATMENT_DATA_END_DATE = 5;

    private static final int RESPONSE_DATA_COUNT = 10;
    private static final int RESPONSE_DATA_PATIENT_ID = 0;
    private static final int RESPONSE_DATA_TIME_POINT = 3;
    private static final int RESPONSE_DATA_DATE = 4;
    private static final int RESPONSE_DATA_RECIST_NOT_DONE = 5;
    private static final int RESPONSE_DATA_RECIST_RESPONSE = 6;
    private static final int RESPONSE_DATA_NO_RECIST_RESPONSE = 7;
    private static final int RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT = 8;
    private static final int RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT_OTHER = 9;

    private WideEcrfFileReader() {
    }

    @NotNull
    public static List<WideFiveDays> readFiveDays(@NotNull String pathToCsv) throws IOException {
        List<WideFiveDays> wideFiveDays = Lists.newArrayList();

        for (CSVRecord line : readCsvSkipHeader(pathToCsv, COMMA_DELIMITER)) {
            if (line.size() == FIVE_DAYS_DATA_COUNT) {
                String widePatientId = WideFileInputInterpreter.toWideID(line.get(FIVE_DAYS_DATA_PATIENT_ID));
                if (!widePatientId.isEmpty()) {
                    boolean dataIsAvailable = line.get(FIVE_DAYS_DATA_IS_AVAILABLE).equals("1");
                    ImmutableWideFiveDays.Builder wideFiveDaysDataBuilder =
                            ImmutableWideFiveDays.builder().widePatientId(widePatientId).dataIsAvailable(dataIsAvailable);

                    if (dataIsAvailable) {
                        wideFiveDaysDataBuilder.informedConsentDate(WideFileInputInterpreter.interpretDateIC(line.get(
                                FIVE_DAYS_DATA_INFORMED_CONSENT_DATE)))
                                .gender(WideFileInputInterpreter.convertGender(line.get(FIVE_DAYS_DATA_GENDER)))
                                .birthYear(WideFileInputInterpreter.interpretBirthYear(line.get(FIVE_DAYS_DATA_BIRTH_YEAR)))
                                .biopsyDate(WideFileInputInterpreter.interpretDateEN(line.get(FIVE_DAYS_DATA_BIOPSY_DATE)))
                                .biopsySite(line.get(FIVE_DAYS_DATA_BIOPSY_SITE))
                                .sampleTissue(line.get(FIVE_DAYS_DATA_SAMPLE_TISSUE))
                                .sampleType(line.get(FIVE_DAYS_DATA_SAMPLE_TYPE))
                                .studyCode(line.get(FIVE_DAYS_DATA_STUDY_CODE))
                                .participatesInOtherTrials(WideFileInputInterpreter.convertParticipatesInOtherTrials(line.get(
                                        FIVE_DAYS_DATA_OTHER_TRIALS)))
                                .otherTrialCodes(line.get(FIVE_DAYS_DATA_OTHER_TRIAL_CODES))
                                .otherTrialStartDates(line.get(FIVE_DAYS_DATA_OTHER_TRIAL_START_DATES));
                    }

                    wideFiveDays.add(wideFiveDaysDataBuilder.build());
                }
            } else {
                LOGGER.warn("Could not properly parse record in WIDE five days csv: {}", line);
            }
        }
        return wideFiveDays;
    }

    @NotNull
    public static List<WidePreAvlTreatmentData> readPreAvlTreatments(@NotNull String pathToCsv) throws IOException {
        List<WidePreAvlTreatmentData> widePreTreatments = Lists.newArrayList();

        for (CSVRecord line : readCsvSkipHeader(pathToCsv, COMMA_DELIMITER)) {
            if (line.size() == PRE_AVL_TREATMENT_DATA_COUNT) {
                String widePatientId = WideFileInputInterpreter.toWideID(line.get(PRE_AVL_TREATMENT_DATA_PATIENT_ID));
                if (!widePatientId.isEmpty()) {
                    WidePreAvlTreatmentData widePreAvlTreatment = ImmutableWidePreAvlTreatmentData.builder()
                            .widePatientId(widePatientId)
                            .hasPreviousTherapy(line.get(PRE_AVL_TREATMENT_DATA_HAS_PREVIOUS_THERAPY).equals("1"))
                            .drug1(line.get(PRE_AVL_TREATMENT_DATA_DRUG1))
                            .drug2(line.get(PRE_AVL_TREATMENT_DATA_DRUG2))
                            .drug3(line.get(PRE_AVL_TREATMENT_DATA_DRUG3))
                            .drug4(line.get(PRE_AVL_TREATMENT_DATA_DRUG4))
                            .lastSystemicTherapyDate(WideFileInputInterpreter.interpretDateEN(line.get(
                                    PRE_AVL_TREATMENT_DATA_LAST_SYSTEMIC_THERAPY_DATE)))
                            .build();

                    widePreTreatments.add(widePreAvlTreatment);
                }
            } else {
                LOGGER.warn("Could not properly parse record in WIDE pre-AVL treatment csv: {}", line);
            }
        }
        return widePreTreatments;
    }

    @NotNull
    public static List<WideBiopsyData> readBiopsies(@NotNull String pathToCsv) throws IOException {
        List<WideBiopsyData> wideBiopsies = Lists.newArrayList();

        for (CSVRecord line : readCsvSkipHeader(pathToCsv, COMMA_DELIMITER)) {
            if (line.size() == BIOPSY_DATA_COUNT) {
                String widePatientId = WideFileInputInterpreter.toWideID(line.get(BIOPSY_DATA_PATIENT_ID));
                if (!widePatientId.isEmpty()) {
                    WideBiopsyData wideBiopsy = ImmutableWideBiopsyData.builder()
                            .widePatientId(widePatientId)
                            .pathologySampleId(line.get(BIOPSY_DATA_PATHOLOGY_SAMPLE_ID))
                            .biopsyDate(WideFileInputInterpreter.interpretDateEN(line.get(BIOPSY_DATA_DATE)))
                            .hasReceivedSuccessfulReport(WideFileInputInterpreter.convertHasReceivedSuccessfulReport(line.get(
                                    BIOPSY_DATA_HAS_RECEIVED_SUCCESSFUL_REPORT)))
                            .build();

                    wideBiopsies.add(wideBiopsy);
                }
            } else {
                LOGGER.warn("Could not properly parse record in WIDE biopsy csv: {}", line);
            }
        }
        return wideBiopsies;
    }

    @NotNull
    private static List<WideAvlTreatmentData> treatmentData(@NotNull List<WideAvlTreatmentData> wideTreatments, CSVRecord line,
            char delimiter) {
        String widePatientId = WideFileInputInterpreter.toWideID(line.get(AVL_TREATMENT_DATA_PATIENT_ID));
        if (!widePatientId.isEmpty()) {

            WideAvlTreatmentData wideTreatment = ImmutableWideAvlTreatmentData.builder()
                    .widePatientId(widePatientId)
                    .drugCode(line.get(AVL_TREATMENT_DATA_DRUG_CODE))
                    .drug(line.get(AVL_TREATMENT_DATA_DRUG))
                    .startDate(delimiter == SEMICOLON_DELIMITER
                            ? WideFileInputInterpreter.interpretDateNL(line.get(AVL_TREATMENT_DATA_START_DATE))
                            : WideFileInputInterpreter.interpretDateEN(line.get(AVL_TREATMENT_DATA_START_DATE)))
                    .endDate(delimiter == SEMICOLON_DELIMITER ? WideFileInputInterpreter.interpretDateNL(line.get(
                            AVL_TREATMENT_DATA_END_DATE)) : WideFileInputInterpreter.interpretDateEN(line.get(AVL_TREATMENT_DATA_END_DATE)))
                    .build();

            wideTreatments.add(wideTreatment);
        }
        return wideTreatments;
    }

    @NotNull
    public static List<WideAvlTreatmentData> readAvlTreatments(@NotNull String pathToCsv) throws IOException {
        List<WideAvlTreatmentData> wideTreatments = Lists.newArrayList();
        char delimiter;
        if (new BufferedReader(new InputStreamReader(new FileInputStream(pathToCsv))).readLine().startsWith("TRAWPTNR,")) {
            delimiter = COMMA_DELIMITER;
        } else if (new BufferedReader(new InputStreamReader(new FileInputStream(pathToCsv))).readLine().startsWith("TRAWPTNR;")) {
            delimiter = SEMICOLON_DELIMITER;
        } else {
            delimiter = ' ';
        }

        for (CSVRecord line : readCsvSkipHeader(pathToCsv, delimiter)) {
            if (line.size() == AVL_TREATMENT_DATA_COUNT) {
                wideTreatments = treatmentData(wideTreatments, line, delimiter);
            } else {
                LOGGER.warn("Could not properly parse record in WIDE treatment csv: {}", line);
            }
        }

        return wideTreatments;
    }

    @NotNull
    public static List<WideResponseData> readResponses(@NotNull String pathToCsv) throws IOException {
        List<WideResponseData> wideResponses = Lists.newArrayList();

        for (CSVRecord line : readCsvSkipHeader(pathToCsv, COMMA_DELIMITER)) {
            if (line.size() == RESPONSE_DATA_COUNT) {
                String widePatientId = WideFileInputInterpreter.toWideID(line.get(RESPONSE_DATA_PATIENT_ID));
                if (!widePatientId.isEmpty()) {
                    boolean recistDone = line.get(RESPONSE_DATA_RECIST_NOT_DONE).equals("FALSE");

                    ImmutableWideResponseData.Builder wideResponseBuilder = ImmutableWideResponseData.builder()
                            .widePatientId(WideFileInputInterpreter.toWideID(line.get(RESPONSE_DATA_PATIENT_ID)))
                            .timePoint(Integer.parseInt(line.get(RESPONSE_DATA_TIME_POINT)))
                            .date(WideFileInputInterpreter.interpretDateEN(line.get(RESPONSE_DATA_DATE)))
                            .recistDone(recistDone)
                            .recistResponse(line.get(RESPONSE_DATA_RECIST_RESPONSE));

                    if (!recistDone) {
                        wideResponseBuilder.noRecistResponse(line.get(RESPONSE_DATA_NO_RECIST_RESPONSE))
                                .noRecistReasonStopTreatment(line.get(RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT))
                                .noRecistReasonStopTreatmentOther(line.get(RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT_OTHER));
                    }

                    wideResponses.add(wideResponseBuilder.build());
                }
            } else {
                LOGGER.warn("Could not properly parse record in WIDE response csv: {}", line);
            }
        }
        return wideResponses;
    }

    @NotNull
    private static List<CSVRecord> readCsvSkipHeader(@NotNull String pathToCsv, char delimiter) throws IOException {
        CSVFormat format = CSVFormat.DEFAULT.withDelimiter(delimiter);
        CSVParser parser = format.parse(new BufferedReader(new InputStreamReader(new FileInputStream(pathToCsv))));

        List<CSVRecord> records = parser.getRecords();
        return records.subList(1, records.size());
    }
}
