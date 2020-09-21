package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class DataUtil {

    private static final String DATE_TIME_FORMAT = "dd-MMM-yyyy";
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    public static final String NONE_STRING = "NONE";
    public static final String NA_STRING = "N/A";

    private DataUtil() {
    }

    @NotNull
    public static String formatPercentage(double percentage) {
        return PERCENTAGE_FORMAT.format(percentage);
    }

    @NotNull
    public static String formatDate(@Nullable LocalDate date) {
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern(DATE_TIME_FORMAT);
        return date != null ? formatter.format(date) : NA_STRING;
    }

    @NotNull
    public static String formatNullableString(@Nullable String string) {
        return string != null ? string : NA_STRING;
    }
}
