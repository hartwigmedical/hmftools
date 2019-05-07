package com.hartwig.hmftools.patientreporter.cfreport;

public final class MathUtil {

    private MathUtil() {
    }

    public static double mapPercentage(final double value, final double inMin, final double inMax) {
        return map(value, inMin, inMax, 0, 100);
    }

    public static double mapPercentageClamped(final double value, final double inMin, final double inMax) {
        return mapClamped(value, inMin, inMax, 0, 100);
    }

    public static double map(final double value, final double inMin, final double inMax, final double outMin, final double outMax) {
        // Remap v from [inMin, inMax] to [outMin, outMax]
        return (value - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
    }

    public static double mapClamped(final double value, final double inMin, final double inMax, final double outMin, final double outMax) {
        return clamp(map(value, inMin, inMax, outMin, outMax), outMin, outMax);
    }

    public static double clamp(final double v, final double min, final double max) {
        return Math.min(max, Math.max(min, v));
    }
}
