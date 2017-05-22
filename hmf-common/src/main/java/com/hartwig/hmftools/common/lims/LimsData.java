package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface LimsData {
    @Nullable
    LocalDate samplingDate();

    @NotNull
    LocalDate arrivalDate();

    boolean isTumor();
}
