package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsBloodData implements LimsData {
    @Nullable
    @Value.Parameter
    public abstract LocalDate samplingDate();

    @NotNull
    @Value.Parameter
    public abstract LocalDate arrivalDate();

    @Override
    public boolean isTumor() {
        return false;
    }
}
