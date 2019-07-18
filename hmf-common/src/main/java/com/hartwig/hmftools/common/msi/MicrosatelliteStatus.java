package com.hartwig.hmftools.common.msi;

import com.hartwig.hmftools.common.numeric.Doubles;

public enum MicrosatelliteStatus {

    MSI,
    MSS,
    UNKNOWN;

    private static final double MSI_CUTOFF = 4.0;

    public static MicrosatelliteStatus fromIndelsPerMb(double microsatelliteIndelsPerMb) {
        return Doubles.greaterOrEqual(microsatelliteIndelsPerMb, MSI_CUTOFF) ? MSI : MSS;
    }

}
