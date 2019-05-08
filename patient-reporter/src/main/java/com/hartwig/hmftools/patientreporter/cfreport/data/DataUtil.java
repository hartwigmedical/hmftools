package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

public final class DataUtil {

    private static final String DATE_TIME_FORMAT = "dd-MMM-yyyy";
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    public static final String NONE_STRING = "NONE";
    public static final String NA_STRING = "N/A";

    private DataUtil() {
    }

    @NotNull
    public static String formatPercentage(final double percentage) {
        return PERCENTAGE_FORMAT.format(percentage);
    }

    @NotNull
    public static String formatDate(@Nullable final LocalDate date) {
        final DateTimeFormatter formatter = DateTimeFormatter.ofPattern(DATE_TIME_FORMAT);
        return date != null ? formatter.format(date) : "?";
    }
}
