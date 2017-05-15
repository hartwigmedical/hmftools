package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class LimsTumorData implements LimsData {
    @Nullable
    private final LocalDate samplingDate;
    @NotNull
    private final LocalDate arrivalDate;
    @Nullable
    private final Double tumorPercentage;

    LimsTumorData(@Nullable final LocalDate samplingDate, @NotNull final LocalDate arrivalDate,
            @Nullable final Double tumorPercentage) {
        this.samplingDate = samplingDate;
        this.arrivalDate = arrivalDate;
        this.tumorPercentage = tumorPercentage;
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

    @Nullable
    public Double tumorPercentage() {
        return tumorPercentage;
    }

    @Override
    public boolean isTumor() {
        return true;
    }
}
