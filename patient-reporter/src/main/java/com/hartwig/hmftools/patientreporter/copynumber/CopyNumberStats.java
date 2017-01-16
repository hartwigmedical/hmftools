package com.hartwig.hmftools.patientreporter.copynumber;

class CopyNumberStats {
    private final int min;
    private final int max;
    private final double mean;

    CopyNumberStats(final int min, final int max, final double mean) {
        this.min = min;
        this.max = max;
        this.mean = mean;
    }

    int min() {
        return min;
    }

    int max() {
        return max;
    }

    double mean() {
        return mean;
    }
}
