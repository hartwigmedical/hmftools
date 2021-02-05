package com.hartwig.hmftools.patientdb.clinical.readers.wide;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.time.format.DateTimeParseException;
import java.util.Locale;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class WideFileInputInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(WideFileInputInterpreter.class);

    private WideFileInputInterpreter() {
    }

    @Nullable
    static LocalDate interpretDateIC(@NotNull String date) {
        return interpretDate(date, DateTimeFormatter.ofPattern("dd/MM/yyyy"));
    }

    @Nullable
    static LocalDate interpretDateEN(@NotNull String date) {
        return interpretDate(date,
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(Locale.ENGLISH));
    }

    @Nullable
    static LocalDate interpretDateNL(@NotNull String date) {
        return interpretDate(date,
                DateTimeFormatter.ofPattern("d-M-yyyy"));
    }

    @Nullable
    private static LocalDate interpretDate(@NotNull String date, @NotNull DateTimeFormatter formatter) {
        if (date.isEmpty()) {
            return null;
        }

        try {
            return LocalDate.parse(date, formatter);
        } catch (DateTimeParseException exception) {
            LOGGER.warn("Could not convert '{}' to date!", date);
            return null;
        }
    }

    @Nullable
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
    static String toWideID(@NotNull String wideIdentifier) {
        return wideIdentifier.replace("-", "");
    }
}
