package com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class FormStatusReader {

    private static final Logger LOGGER = LogManager.getLogger(FormStatusReader.class);

    private static final int PATIENT_KEY_COLUMN = 0;
    private static final int PATIENT_ID_COLUMN = 1;
    private static final int DEPT_PERS_ID_COLUMN = 2;
    private static final int EMAIL_COLUMN = 3;
    private static final int INSTITUTE_ID_COLUMN = 4;
    private static final int INSTITUTE_NAME_COLUMN = 5;
    private static final int COMPLETED_FORM_ID_COLUMN = 6;
    private static final int DATA_STATUS_COLUMN = 7;
    private static final int FORM_ID_COLUMN = 8;
    private static final int FORM_NAME_COLUMN = 9;
    private static final int FORM_SEQ_NUM_COLUMN = 10;
    private static final int DATA_STATUS_STRING_COLUMN = 11;
    private static final int STUDY_EVENT_ID_COLUMN = 12;
    private static final int STUDY_EVENT_SEQ_NUM_COLUMN = 13;
    private static final int STUDY_EVENT_NAME_COLUMN = 14;
    private static final int LAST_SAVED_COLUMN = 15;
    private static final int LOCKED_COLUMN = 16;

    private static final int FIELD_COUNT = 17;

    private static final String FIELD_SEPARATOR = ",";
    private static final char QUOTE = '"';

    private FormStatusReader() {
    }

    @NotNull
    public static FormStatusModel buildModelFromCsv(@NotNull final String pathToCsv) throws IOException {
        Map<FormStatusKey, FormStatus> formStatuses = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(pathToCsv).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = splitCsvLine(line, FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                FormStatusKey formKey = new ImmutableFormStatusKey(removeQuotes(parts[PATIENT_ID_COLUMN].replaceAll("-", "")),
                        removeParentheses(removeQuotes(parts[FORM_NAME_COLUMN])),
                        removeQuotes(parts[FORM_SEQ_NUM_COLUMN]),
                        removeParentheses(removeQuotes(parts[STUDY_EVENT_NAME_COLUMN])),
                        removeQuotes(parts[STUDY_EVENT_SEQ_NUM_COLUMN]));
                FormStatus formStatus = new ImmutableFormStatus(interpretState(removeQuotes(parts[DATA_STATUS_COLUMN])),
                        interpretLocked(removeQuotes(parts[LOCKED_COLUMN])));
                formStatuses.put(formKey, formStatus);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in form status csv: {}", line);
            }
        }
        return new ImmutableFormStatusModel(formStatuses);
    }

    @NotNull
    private static String[] splitCsvLine(@NotNull String line, @NotNull String separator, int limit) {
        // Split a line in a csv/tsv file by the separator, ignoring occurrences of separator in quoted groups
        return line.split(separator + "(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", limit);
    }

    @NotNull
    private static String removeQuotes(@NotNull String text) {
        String trimmedText = text.trim();
        if (trimmedText.length() > 0 && trimmedText.charAt(0) == QUOTE && trimmedText.charAt(trimmedText.length() - 1) == QUOTE) {
            return trimmedText.substring(1, trimmedText.length() - 1);
        }
        return trimmedText;
    }

    @NotNull
    private static String removeParentheses(@NotNull String text) {
        final Pattern pattern = Pattern.compile(".*(?=(?:\\([0-9]+\\))$)");
        final Matcher matcher = pattern.matcher(text);
        if (matcher.find()) {
            return matcher.group(0).trim();
        } else {
            return text;
        }
    }

    private static boolean interpretLocked(@NotNull String locked) {
        if (!locked.equalsIgnoreCase("true") && !locked.equalsIgnoreCase("false")) {
            LOGGER.warn("Undefined locked status - cannot convert to boolean: {}", locked);
        }

        return locked.equalsIgnoreCase("true");
    }

    @NotNull
    private static FormStatusState interpretState(@NotNull String state) {
        switch (state) {
            case "0":
                return FormStatusState.SAVED;
            case "1":
                return FormStatusState.SUBMITTED;
            case "2":
                return FormStatusState.SUBMITTED_WITH_DISCREPANCIES;
            case "3":
                return FormStatusState.SUBMITTED_WITH_MISSING;
            case "4":
                return FormStatusState.VERIFIED;
            default: {
                LOGGER.warn("Could not interpret data status: {}", state);
                return FormStatusState.UNDEFINED;
            }
        }
    }
}
