package com.hartwig.hmftools.patientdb.dao;

import org.jetbrains.annotations.Nullable;

final public class DatabaseUtil {

    @Nullable
    public static Double decimal(@Nullable Double number) {
        return number == null ? null : decimal(number.doubleValue());
    }

    public static double decimal(double number) {
        return round(number, 4);
    }

    private static double round(double number, int decimalPoints) {
        double multiplier = Math.pow(10, decimalPoints);
        return Math.round(number * multiplier) / multiplier;
    }
}
