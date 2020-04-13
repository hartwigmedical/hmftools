package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeFormatterBuilder;
import java.util.Locale;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class WideFileInputInterpreter {

    private WideFileInputInterpreter() {
    }

    @Nullable
    public static LocalDate interpretDateIC(@NotNull String date) {
        if (date.isEmpty()) {
            return null;
        }

        return LocalDate.parse(date, DateTimeFormatter.ofPattern("dd/MM/yyyy"));
    }

    @Nullable
    public static LocalDate interpretDateEN(@NotNull String date) {
        if (date.isEmpty()) {
            return null;
        }

        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(Locale.ENGLISH);
        return LocalDate.parse(date, inputFormatter);
    }

    @Nullable
    public static LocalDate interpretDateNL(@NotNull String date) {
        if (date.isEmpty()) {
            return null;
        }

        DateTimeFormatter inputFormatter =
                new DateTimeFormatterBuilder().parseCaseInsensitive().appendPattern("dd-MMM-yyyy").toFormatter(new Locale("nl", "NL"));
        return LocalDate.parse(date, inputFormatter);
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
