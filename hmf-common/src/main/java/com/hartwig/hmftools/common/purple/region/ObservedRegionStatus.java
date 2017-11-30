package com.hartwig.hmftools.common.purple.region;

import org.jetbrains.annotations.NotNull;

public enum ObservedRegionStatus {
    GERMLINE_HOM_DELETION,
    GERMLINE_HET_DELETION,
    GERMLINE_AMPLIFICATION,
    GERMLINE_NOISE,
    CLUSTER,
    DIPLOID,
    UNKNOWN;

    public static ObservedRegionStatus fromString(@NotNull final String status) {
        if (status.equals("SOMATIC")) {
            return DIPLOID;
        }

        return ObservedRegionStatus.valueOf(status);
    }
}
