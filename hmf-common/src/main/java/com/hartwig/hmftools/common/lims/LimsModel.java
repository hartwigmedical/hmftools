package com.hartwig.hmftools.common.lims;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LimsModel {

    @NotNull
    private final Map<String, LimsData> dataPerSample;

    LimsModel(@NotNull final Map<String, LimsData> dataPerSample) {
        this.dataPerSample = dataPerSample;
    }

    @NotNull
    @VisibleForTesting
    Map<String, LimsData> data() {
        return dataPerSample;
    }

    @Nullable
    public LimsData findDataPerSample(@NotNull final String sample) {
        return dataPerSample.get(sample);
    }
}
