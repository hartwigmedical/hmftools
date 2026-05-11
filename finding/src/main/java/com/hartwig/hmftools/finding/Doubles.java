package com.hartwig.hmftools.finding;

import java.math.BigDecimal;
import java.math.RoundingMode;

@SuppressWarnings("unused")
final class Doubles
{
    private static final double EPSILON = 1e-10;

    public static boolean equal(double first, double second)
    {
        return Math.abs(first - second) < EPSILON;
    }

    public static boolean lessThan(double value, double reference)
    {
        return value - reference < -EPSILON;
    }

    public static boolean lessOrEqual(double value, double reference)
    {
        return value - reference < EPSILON;
    }

    public static boolean greaterThan(double value, double reference)
    {
        return value - reference > EPSILON;
    }

    public static boolean greaterOrEqual(double value, double reference)
    {
        return value - reference > -EPSILON;
    }

    public static double round(double value, int decimals) {
        return toBigDecimalRounded(value, decimals).doubleValue();
    }

    private static BigDecimal toBigDecimalRounded(double value, int decimals) {
        return new BigDecimal(value).setScale(decimals, RoundingMode.HALF_UP);
    }
}
