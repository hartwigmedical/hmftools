package com.hartwig.hmftools.patientdb.readers.wide;

import java.io.File;
import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.List;
import java.util.Locale;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class WideEcrfFileReader {

    private static final Logger LOGGER = LogManager.getLogger(WideEcrfFileReader.class);
    private static final String SEMICOLON = ";";
    private static final String COMMA = ",";

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

    private static final int PRE_AVL_TREATMENT_DATA_COUNT = 8;
    private static final int PRE_AVL_TREATMENT_DATA_PATIENT_ID = 0;
    private static final int PRE_AVL_TREATMENT_DATA_HAS_PREVIOUS_THERAPY = 1;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG1 = 2;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG2 = 3;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG3 = 4;
    private static final int PRE_AVL_TREATMENT_DATA_DRUG4 = 5;
    private static final int PRE_AVL_TREATMENT_DATA_LAST_SYSTEMIC_THERAPY_DATE = 6;

    private static final int BIOPSY_DATA_COUNT = 6;
    private static final int BIOPSY_DATA_PATIENT_ID = 1;
    private static final int BIOPSY_DATA_PATHOLOGY_SAMPLE_ID = 3;
    private static final int BIOPSY_DATA_DATE = 4;
    private static final int BIOPSY_DATA_HAS_RECEIVED_SUCCESSFUL_REPORT = 5;

    private static final int AVL_TREATMENT_DATA_COUNT = 5;
    private static final int AVL_TREATMENT_DATA_PATIENT_ID = 0;
    private static final int AVL_TREATMENT_DATA_DRUG_CODE = 1;
    private static final int AVL_TREATMENT_DATA_DRUG = 2;
    private static final int AVL_TREATMENT_DATA_START_DATE = 3;
    private static final int AVL_TREATMENT_DATA_END_DATE = 4;

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

        List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(COMMA, FIVE_DAYS_DATA_COUNT);
            if (parts.length == FIVE_DAYS_DATA_COUNT) {
                WideFiveDays wideFiveDaysData = ImmutableWideFiveDays.builder()
                        .patientId(toWideID(parts[FIVE_DAYS_DATA_PATIENT_ID]))
                        .dataIsAvailable(parts[FIVE_DAYS_DATA_IS_AVAILABLE].equals("1"))
                        .informedConsentDate(interpretDateIC(parts[FIVE_DAYS_DATA_INFORMED_CONSENT_DATE]))
                        .gender(convertGender(parts[FIVE_DAYS_DATA_GENDER]))
                        .birthYear(interpretBirthYear(parts[FIVE_DAYS_DATA_BIRTH_YEAR]))
                        .biopsyDate(interpretDateEN(parts[FIVE_DAYS_DATA_BIOPSY_DATE]))
                        .biopsySite(parts[FIVE_DAYS_DATA_BIOPSY_SITE])
                        .sampleTissue(parts[FIVE_DAYS_DATA_SAMPLE_TISSUE])
                        .sampleType(parts[FIVE_DAYS_DATA_SAMPLE_TYPE])
                        .studyCode(parts[FIVE_DAYS_DATA_STUDY_CODE])
                        .participatesInOtherTrials(convertParticipatesInOtherTrials(parts[FIVE_DAYS_DATA_OTHER_TRIALS]))
                        .otherTrialCodes(parts[FIVE_DAYS_DATA_OTHER_TRIAL_CODES])
                        .otherTrialStartDates(parts[FIVE_DAYS_DATA_OTHER_TRIAL_START_DATES])
                        .build();

                wideFiveDays.add(wideFiveDaysData);
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
            String[] parts = line.split(COMMA, PRE_AVL_TREATMENT_DATA_COUNT);
            if (parts.length == PRE_AVL_TREATMENT_DATA_COUNT) {
                WidePreAvlTreatmentData widePreAvlTreatment = ImmutableWidePreAvlTreatmentData.builder()
                        .patientId(toWideID(parts[PRE_AVL_TREATMENT_DATA_PATIENT_ID]))
                        .hasPreviousTherapy(parts[PRE_AVL_TREATMENT_DATA_HAS_PREVIOUS_THERAPY].equals("1"))
                        .drug1(parts[PRE_AVL_TREATMENT_DATA_DRUG1])
                        .drug2(parts[PRE_AVL_TREATMENT_DATA_DRUG2])
                        .drug3(parts[PRE_AVL_TREATMENT_DATA_DRUG3])
                        .drug4(parts[PRE_AVL_TREATMENT_DATA_DRUG4])
                        .lastSystemicTherapyDate(interpretDateEN(parts[PRE_AVL_TREATMENT_DATA_LAST_SYSTEMIC_THERAPY_DATE]))
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
            String[] parts = line.split(COMMA, BIOPSY_DATA_COUNT);
            if (parts.length == BIOPSY_DATA_COUNT) {
                WideBiopsyData wideBiopsy = ImmutableWideBiopsyData.builder()
                        .patientId(toWideID(parts[BIOPSY_DATA_PATIENT_ID]))
                        .pathologySampleId(parts[BIOPSY_DATA_PATHOLOGY_SAMPLE_ID])
                        .biopsyDate(interpretDateEN(parts[BIOPSY_DATA_DATE]))
                        .hasReceivedSuccessfulReport(convertHasReceivedSuccessfulReport(parts[BIOPSY_DATA_HAS_RECEIVED_SUCCESSFUL_REPORT]))
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
            String[] parts = line.split(SEMICOLON, AVL_TREATMENT_DATA_COUNT);
            if (parts.length == AVL_TREATMENT_DATA_COUNT) {
                WideAvlTreatmentData wideTreatment = ImmutableWideAvlTreatmentData.builder()
                        .patientId(toWideID(parts[AVL_TREATMENT_DATA_PATIENT_ID]))
                        .drugCode(parts[AVL_TREATMENT_DATA_DRUG_CODE])
                        .drug(parts[AVL_TREATMENT_DATA_DRUG])
                        .startDate(interpretDateEN(parts[AVL_TREATMENT_DATA_START_DATE]))
                        .endDate(interpretDateEN(parts[AVL_TREATMENT_DATA_END_DATE]))
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
            String[] parts = line.split(COMMA, RESPONSE_DATA_COUNT);
            if (parts.length == RESPONSE_DATA_COUNT) {
                WideResponseData wideResponse = ImmutableWideResponseData.builder()
                        .patientId(toWideID(parts[RESPONSE_DATA_PATIENT_ID]))
                        .timePoint(Integer.parseInt(parts[RESPONSE_DATA_TIME_POINT]))
                        .date(interpretDateNL(parts[RESPONSE_DATA_DATE]))
                        .recistDone(parts[RESPONSE_DATA_RECIST_NOT_DONE].equals("FALSE"))
                        .recistResponse(parts[RESPONSE_DATA_RECIST_RESPONSE])
                        .noRecistResponse(parts[RESPONSE_DATA_NO_RECIST_RESPONSE])
                        .noRecistReasonStopTreatment(parts[RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT])
                        .noRecistReasonStopTreatmentOther(parts[RESPONSE_DATA_NO_RECIST_REASON_STOP_TREATMENT_OTHER])
                        .build();

                wideResponses.add(wideResponse);
            } else {
                LOGGER.warn("Could not properly parse line in WIDE response csv: {}", line);
            }
        }
        return wideResponses;
    }

    @Nullable
    @VisibleForTesting
    public static LocalDate interpretDateIC(@NotNull String date) {
        if (date.isEmpty()) {
            return null;
        }

        return LocalDate.parse(date, DateTimeFormatter.ofPattern("dd/MM/yyyy"));
    }

    @Nullable
    @VisibleForTesting
    public static LocalDate interpretDateEN(@NotNull String date) {
        if (date.isEmpty()) {
            return null;
        }

        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(Locale.ENGLISH);
        LocalDate localDate = LocalDate.parse(date, inputFormatter);

        DateTimeFormatter outputFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
        String formattedString = localDate.format(outputFormatter);
        return LocalDate.parse(formattedString);
    }

    @Nullable
    @VisibleForTesting
    public static LocalDate interpretDateNL(@NotNull String date) {
        if (date.isEmpty()) {
            return null;
        }

        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(new Locale("nl", "NL"));
        LocalDate localDate = LocalDate.parse(date, inputFormatter);

        DateTimeFormatter outputFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");
        String formattedString = localDate.format(outputFormatter);
        return LocalDate.parse(formattedString);
    }

    @Nullable
    @VisibleForTesting
    static Integer interpretBirthYear(@NotNull String birthYear) {
        if (birthYear.isEmpty()) {
            return null;
        }

        return Integer.parseInt(birthYear);
    }

    @Nullable
    @VisibleForTesting
    static String convertGender(@NotNull String gender) {
        if (gender.equals("1")) {
            return "male";
        } else if (gender.equals("2")) {
            return "female";
        } else {
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    static Boolean convertParticipatesInOtherTrials(@NotNull String participatesInOtherTrials) {
        if (participatesInOtherTrials.equals("Y")) {
            return true;
        } else if (participatesInOtherTrials.equals("N")) {
            return false;
        } else {
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    static Boolean convertHasReceivedSuccessfulReport(@NotNull String hasReceivedSuccessfulReport) {
        if (hasReceivedSuccessfulReport.equals("yes")) {
            return true;
        } else if (hasReceivedSuccessfulReport.equals("no")) {
            return false;
        } else {
            return null;
        }
    }

    @NotNull
    @VisibleForTesting
    static String toWideID(@NotNull String wideIdentifier) {
        return wideIdentifier.replace("-", "");
    }
}
