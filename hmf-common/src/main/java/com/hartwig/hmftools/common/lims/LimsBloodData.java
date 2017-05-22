package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class LimsBloodData implements LimsData {
    @Nullable
    private final LocalDate samplingDate;
    @NotNull
    private final LocalDate arrivalDate;

    LimsBloodData(@Nullable final LocalDate samplingDate, @NotNull final LocalDate arrivalDate) {
        this.samplingDate = samplingDate;
        this.arrivalDate = arrivalDate;
    }

    @Nullable
    @Override
    public LocalDate samplingDate() {
        return samplingDate;
    }

    @NotNull
    @Override
    public LocalDate arrivalDate() {
        return arrivalDate;
    }

    @Override
    public boolean isTumor() {
        return false;
    }
}
