package com.hartwig.hmftools.patientdb.dao;

final public class DatabaseUtil {

    public static double decimal(double number) {
        return round(number, 4);
    }

    private static double round(double number, int decimalPoints) {
        double multiplier = Math.pow(10, decimalPoints);
        return Math.round(number * multiplier) / multiplier;
    }
}
