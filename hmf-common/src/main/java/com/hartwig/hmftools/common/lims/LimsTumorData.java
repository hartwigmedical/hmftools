package com.hartwig.hmftools.common.lims;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class LimsTumorData implements LimsData {
    @Nullable
    @Value.Parameter
    public abstract LocalDate samplingDate();

    @NotNull
    @Value.Parameter
    public abstract LocalDate arrivalDate();

    @Nullable
    @Value.Parameter
    public abstract Double tumorPercentage();

    @Override
    public boolean isTumor() {
        return true;
    }
}
