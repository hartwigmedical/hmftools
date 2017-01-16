package com.hartwig.hmftools.patientreporter.copynumber;

public class CopyNumberStats {
    private final int min;
    private final int max;
    private final double mean;

    CopyNumberStats(final int min, final int max, final double mean) {
        this.min = min;
        this.max = max;
        this.mean = mean;
    }

    public int min() {
        return min;
    }

    public int max() {
        return max;
    }

    public double mean() {
        return mean;
    }
}
