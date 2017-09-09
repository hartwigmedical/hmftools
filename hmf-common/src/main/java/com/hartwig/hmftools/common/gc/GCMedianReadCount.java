package com.hartwig.hmftools.common.gc;

import org.jetbrains.annotations.NotNull;

public interface GCMedianReadCount {

    int meanReadCount();

    int medianReadCount();

    int medianReadCount(@NotNull GCBucket bucket);

    default int medianReadCount(@NotNull GCProfile profile) {
        return medianReadCount(GCBucket.create(profile));
    }
}
