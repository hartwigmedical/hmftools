package com.hartwig.hmftools.patientdb.readers.wide;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class WideEcrfFileReader {

    private static final Logger LOGGER = LogManager.getLogger(WideEcrfFileReader.class);
    private static final String SEMICOLON = ";";
    private static final String COMMA = ",";

    private static final int FIVE_DAYS_DATA_MIN_COUNT = 7;
    private static final int FIVE_DAYS_DATA_MAX_COUNT = 15;
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

    private static final int BIOPSY_DATA_MIN_COUNT = 5;
    private static final int BIOPSY_DATA_MAX_COUNT = 6;
    private static final int BIOPSY_DATA_PATIENT_ID = 1;
    private static final int BIOPSY_DATA_PATHOLOGY_SAMPLE_ID = 3;
    private static final int BIOPSY_DATA_DATE = 4;
    private static final int BIOPSY_DATA_HAS_RECEIVED_SUCCESSFUL_REPORT = 5;

    private static final int AVL_TREATMENT_DATA_MIN_COUNT = 5;
    private static final int AVL_TREATMENT_DATA_MAX_COUNT = 6;
    private static final int AVL_TREATMENT_DATA_PATIENT_ID = 0;
    private static final int AVL_TREATMENT_DATA_DRUG_CODE = 2;
    private static final int AVL_TREATMENT_DATA_DRUG = 3;
    private static final int AVL_TREATMENT_DATA_START_DATE = 4;
    private static final int AVL_TREATMENT_DATA_END_DATE = 5;

    private static final int RESPONSE_DATA_MIN_COUNT = 6;
    private static final int RESPONSE_DATA_MAX_COUNT = 10;
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

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(COMMA);
            if (parts.length >= FIVE_DAYS_DATA_MIN_COUNT && parts.length <= FIVE_DAYS_DATA_MAX_COUNT) {
                boolean dataIsAvailable = parts[FIVE_DAYS_DATA_IS_AVAILABLE].equals("1");
                ImmutableWideFiveDays.Builder wideFiveDaysDataBuilder = ImmutableWideFiveDays.builder()
                        .patientId(WideFileInputInterpreter.toWideID(parts[FIVE_DAYS_DATA_PATIENT_ID]))
                        .dataIsAvailable(dataIsAvailable);

                if (dataIsAvailable) {
                    wideFiveDaysDataBuilder.informedConsentDate(WideFileInputInterpreter.interpretDateIC(parts[FIVE_DAYS_DATA_INFORMED_CONSENT_DATE]))
                            .gender(WideFileInputInterpreter.convertGender(parts[FIVE_DAYS_DATA_GENDER]))
                            .birthYear(WideFileInputInterpreter.interpretBirthYear(parts[FIVE_DAYS_DATA_BIRTH_YEAR]))
                            .biopsyDate(WideFileInputInterpreter.interpretDateEN(parts[FIVE_DAYS_DATA_BIOPSY_DATE]))
                            .biopsySite(parts[FIVE_DAYS_DATA_BIOPSY_SITE])
                            .sampleTissue(parts[FIVE_DAYS_DATA_SAMPLE_TISSUE])
                            .sampleType(parts[FIVE_DAYS_DATA_SAMPLE_TYPE])
                            .studyCode(parts[FIVE_DAYS_DATA_STUDY_CODE])
                            .participatesInOtherTrials(WideFileInputInterpreter.convertParticipatesInOtherTrials(parts[FIVE_DAYS_DATA_OTHER_TRIALS]))
                            .otherTrialCodes(parts.length > FIVE_DAYS_DATA_OTHER_TRIAL_CODES
                                    ? parts[FIVE_DAYS_DATA_OTHER_TRIAL_CODES]
                                    : Strings.EMPTY)
                            .otherTrialStartDates(parts.length > FIVE_DAYS_DATA_OTHER_TRIAL_START_DATES
                                    ? parts[FIVE_DAYS_DATA_OTHER_TRIAL_START_DATES]
                                    : Strings.EMPTY);
                }

                wideFiveDays.add(wideFiveDaysDataBuilder.build());
            } else {
                LOGGER.warn("Could not properly parse line in WIDE five days csv: {}", line);
            }
        }
        return wideFiveDays;
    }

    @NotNull
    public static List<WidePreAvlTreatmentData> readPreAvlTreatments(@NotNull String pathToCsv) throws IOException {
        List<WidePreAvlTreatmentData> widePreTreatments = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(COMMA);
            if (parts.length == PRE_AVL_TREATMENT_DATA_COUNT) {
                WidePreAvlTreatmentData widePreAvlTreatment = ImmutableWidePreAvlTreatmentData.builder()
                        .patientId(WideFileInputInterpreter.toWideID(parts[PRE_AVL_TREATMENT_DATA_PATIENT_ID]))
                        .hasPreviousTherapy(parts[PRE_AVL_TREATMENT_DATA_HAS_PREVIOUS_THERAPY].equals("1"))
                        .drug1(parts[PRE_AVL_TREATMENT_DATA_DRUG1])
                        .drug2(parts[PRE_AVL_TREATMENT_DATA_DRUG2])
                        .drug3(parts[PRE_AVL_TREATMENT_DATA_DRUG3])
                        .drug4(parts[PRE_AVL_TREATMENT_DATA_DRUG4])
                        .lastSystemicTherapyDate(WideFileInputInterpreter.interpretDateEN(parts[PRE_AVL_TREATMENT_DATA_LAST_SYSTEMIC_THERAPY_DATE]))
                        .build();

                widePreTreatments.add(widePreAvlTreatment);
            } else {
                LOGGER.warn("Could not properly parse line in WIDE pre-AVL treatment csv: {}", line);
            }
        }
        return widePreTreatments;
    }

    @NotNull
    public static List<WideBiopsyData> readBiopsies(@NotNull String pathToCsv) throws IOException {
        List<WideBiopsyData> wideBiopsies = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(COMMA);
            if (parts.length >= BIOPSY_DATA_MIN_COUNT && parts.length <= BIOPSY_DATA_MAX_COUNT) {
                Boolean hasSuccessfulReport = parts.length > BIOPSY_DATA_HAS_RECEIVED_SUCCESSFUL_REPORT ?
                        WideFileInputInterpreter.convertHasReceivedSuccessfulReport(parts[BIOPSY_DATA_HAS_RECEIVED_SUCCESSFUL_REPORT]) :
                        null;
                WideBiopsyData wideBiopsy = ImmutableWideBiopsyData.builder()
                        .patientId(WideFileInputInterpreter.toWideID(parts[BIOPSY_DATA_PATIENT_ID]))
                        .pathologySampleId(parts[BIOPSY_DATA_PATHOLOGY_SAMPLE_ID])
                        .biopsyDate(WideFileInputInterpreter.interpretDateEN(parts[BIOPSY_DATA_DATE]))
                        .hasReceivedSuccessfulReport(hasSuccessfulReport)
                        .build();

                wideBiopsies.add(wideBiopsy);
            } else {
                LOGGER.warn("Could not properly parse line in WIDE biopsy csv: {}", line);
            }
        }
        return wideBiopsies;
    }

    @NotNull
    public static List<WideAvlTreatmentData> readAvlTreatments(@NotNull String pathToCsv) throws IOException {
        List<WideAvlTreatmentData> wideTreatments = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(SEMICOLON);
            if (parts.length >= AVL_TREATMENT_DATA_MIN_COUNT && parts.length <= AVL_TREATMENT_DATA_MAX_COUNT) {
                // End date is optional
                LocalDate endDate = parts.length > AVL_TREATMENT_DATA_END_DATE
                        ? WideFileInputInterpreter.interpretDateEN(parts[AVL_TREATMENT_DATA_END_DATE])
                        : null;
                WideAvlTreatmentData wideTreatment = ImmutableWideAvlTreatmentData.builder()
                        .patientId(WideFileInputInterpreter.toWideID(parts[AVL_TREATMENT_DATA_PATIENT_ID]))
                        .drugCode(parts[AVL_TREATMENT_DATA_DRUG_CODE])
                        .drug(parts[AVL_TREATMENT_DATA_DRUG])
                        .startDate(WideFileInputInterpreter.interpretDateEN(parts[AVL_TREATMENT_DATA_START_DATE]))
                        .endDate(endDate)
                        .build();

                wideTreatments.add(wideTreatment);
            } else {
                LOGGER.warn("Could not properly parse line in WIDE AVL treatment csv: {}", line);
            }
        }
        return wideTreatments;
    }

    @NotNull
    public static List<WideResponseData> readResponses(@NotNull String pathToCsv) throws IOException {
        List<WideResponseData> wideResponses = Lists.newArrayList();

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(COMMA);
            if (parts.length >= RESPONSE_DATA_MIN_COUNT && parts.length <= RESPONSE_DATA_MAX_COUNT) {
                boolean recistDone = parts[RESPONSE_DATA_RECIST_NOT_DONE].equals("FALSE");
                ImmutableWideResponseData.Builder wideResponseBuilder = ImmutableWideResponseData.builder()
                        .patientId(WideFileInputInterpreter.toWideID(parts[RESPONSE_DATA_PATIENT_ID]))
                        .timePoint(Integer.parseInt(parts[RESPONSE_DATA_TIME_POINT]))
                        .date(WideFileInputInterpreter.interpretDateNL(parts[RESPONSE_DATA_DATE]))
                        .recistDone(recistDone)
                        .recistResponse(
                                parts.length > RESPONSE_DATA_RECIST_RESPONSE ? parts[RESPONSE_DATA_RECIST_RESPONSE] : Strings.EMPTY);

                if (!recistDone) {
                    wideResponseBuilder.noRecistResponse(
                            parts.length > RESPONSE_DATA_NO_RECIST_RESPONSE ? parts[RESPONSE_DATA_NO_RECIST_RESPONSE] : Strings.EMPTY)
                            .noRecistReasonStopTreatment(parts.length > RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT
                                    ? parts[RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT]
                                    : Strings.EMPTY)
                            .noRecistReasonStopTreatmentOther(parts.length > RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT_OTHER
                                    ? parts[RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT_OTHER]
                                    : Strings.EMPTY);
                }

                wideResponses.add(wideResponseBuilder.build());
            } else {
                LOGGER.warn("Could not properly parse line in WIDE response csv: {}", line);
            }
        }
        return wideResponses;
    }
}
