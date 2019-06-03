package com.hartwig.hmftools.patientdb.dao;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.TableField;

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

    @NotNull
    public static String checkStringLength(@NotNull String str, @NotNull TableField<?, String> field) {
        int maxLength = field.getDataType().length();
        return str.length() > maxLength ? str.substring(0, maxLength) : str;
    }

    public static double getValueNotNull(@Nullable Double value) {
        return value != null ? value : 0D;
    }

    public static int getValueNotNull(@Nullable Integer value) {
        return value != null ? value : 0;
    }

    public static byte getValueNotNull(@Nullable Byte value) {
        return value != null ? value : 0;
    }

    @NotNull
    public static String getValueNotNull(@Nullable String value) {
        return value != null ? value : Strings.EMPTY;
    }

}
