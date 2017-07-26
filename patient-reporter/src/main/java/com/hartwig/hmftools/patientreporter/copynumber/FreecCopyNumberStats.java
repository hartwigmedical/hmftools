package com.hartwig.hmftools.patientreporter.copynumber;

class FreecCopyNumberStats {
    private final int min;
    private final int max;
    private final double mean;

    FreecCopyNumberStats(final int min, final int max, final double mean) {
        this.min = min;
        this.max = max;
        this.mean = mean;
    }

    int min() {
        return min;
    }

    @SuppressWarnings("unused")
    int max() {
        return max;
    }

    @SuppressWarnings("unused")
    double mean() {
        return mean;
    }
}
