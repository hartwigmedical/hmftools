package com.hartwig.hmftools.patientreporter.cfreport;

public final class MathUtil {

    private MathUtil() {
    }

    public static double mapPercentage(double value, double inMin, double inMax) {
        return map(value, inMin, inMax, 0, 100);
    }

    public static double map(double value, double inMin, double inMax, double outMin, double outMax) {
        // Remap v from [inMin, inMax] to [outMin, outMax]
        return (value - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
    }

    public static double mapClamped(double value, double inMin, double inMax, double outMin, double outMax) {
        return clamp(map(value, inMin, inMax, outMin, outMax), outMin, outMax);
    }

    public static double clamp(double v, double min, double max) {
        return Math.min(max, Math.max(min, v));
    }
}
