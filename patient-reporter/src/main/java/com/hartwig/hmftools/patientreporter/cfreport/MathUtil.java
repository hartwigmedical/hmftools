package com.hartwig.hmftools.patientreporter.cfreport;

public final class MathUtil {

    /**
     * Remap v from [inMin, inMax] to [0, 100]
     */
    public static double mapPercentage(final double v, final double inMin, final double inMax) {
        return map(v, inMin, inMax, 0, 100);
    }

    /**
     * Remap v from [inMin, inMax] to [0, 100]
     */
    public static double mapPercentageClamped(final double v, final double inMin, final double inMax) {
        return mapClamped(v, inMin, inMax, 0, 100);
    }

    /**
     * Remap v from [inMin, inMax] to [outMin, outMax]
     */
    public static double map(final double v, final double inMin, final double inMax, final double outMin, final double outMax) {
        return (v - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
    }

    /**
     * Remap v from [inMin, inMax] to [outMin, outMax]
     */
    public static double mapClamped(final double v, final double inMin, final double inMax, final double outMin, final double outMax) {
        return clamp(map(v, inMin, inMax, outMin, outMax), outMin, outMax);
    }

    /**
     * Limit v to [min, max]
     */
    public static double clamp(final double v, final double min, final double max) {
        return Math.min(max, Math.max(min, v));
    }

}
