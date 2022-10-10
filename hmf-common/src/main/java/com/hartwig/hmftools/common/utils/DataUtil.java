package com.hartwig.hmftools.common.utils;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class DataUtil {

    private static final DateTimeFormatter DATE_TIME_FORMATTER = DateTimeFormatter.ofPattern("dd-MMM-yyyy");
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");
    private static final DecimalFormat PERCENTAGE_FORMAT_WITH_DIGIT = new DecimalFormat("#.#'%'");
    public static final String NONE_STRING = "NONE";
    public static final String NA_STRING = "N/A";

    private DataUtil() {
    }

    @NotNull
    public static String formatPercentageDigit(double percentage) {
        return PERCENTAGE_FORMAT_WITH_DIGIT.format(percentage * 100);
    }

    @NotNull
    public static String formatPercentage(double percentage) {
        return PERCENTAGE_FORMAT.format(percentage);
    }

    @NotNull
    public static String formatDate(@Nullable LocalDate date) {
        return date != null ? DATE_TIME_FORMATTER.format(date) : NA_STRING;
    }

    @NotNull
    public static String formatNullableString(@Nullable String string) {
        return string != null ? string : NA_STRING;
    }
}
