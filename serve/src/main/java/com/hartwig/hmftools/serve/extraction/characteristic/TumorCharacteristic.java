package com.hartwig.hmftools.serve.extraction.characteristic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorCharacteristic {

    @Nullable
    public abstract TumorCharacteristicAnnotation tumorCharacteristicAnnotation();

    @Nullable
    public abstract TumorCharacteristicsAtLeast atLeast();

    @Nullable
    public abstract Double cutOff();
}
