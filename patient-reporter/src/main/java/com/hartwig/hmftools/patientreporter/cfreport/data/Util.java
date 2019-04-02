package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.text.DecimalFormat;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;

public class Util {

   // Number formatting
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    // Missing/invalid data indicators
    public static final String NoneString = "NONE";
    public static final String NAString = "N/A";

    @NotNull
    public static String formatPercentage(final double percentage) {
        return PERCENTAGE_FORMAT.format(percentage);
    }

    @NotNull
    public static String formatDate(@Nullable final LocalDate date) {
        final DateTimeFormatter formatter = DateTimeFormatter.ofPattern(ReportResources.DATE_TIME_FORMAT);
        return date != null ? formatter.format(date) : "?";
    }

}
