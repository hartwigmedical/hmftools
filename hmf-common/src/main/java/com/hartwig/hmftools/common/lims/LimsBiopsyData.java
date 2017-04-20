package com.hartwig.hmftools.common.lims;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class LimsBiopsyData {

    @Nullable
    private final String samplingDate;
    @NotNull
    private final String arrivalDate;
    private final double tumorPercentage;

    LimsBiopsyData(@Nullable final String samplingDate, @NotNull final String arrivalDate,
            final double tumorPercentage) {
        this.samplingDate = samplingDate;
        this.arrivalDate = arrivalDate;
        this.tumorPercentage = tumorPercentage;
    }

    @Nullable
    String samplingDate() {
        return samplingDate;
    }

    @NotNull
    String arrivalDate() {
        return arrivalDate;
    }

    double tumorPercentage() {
        return tumorPercentage;
    }
}
