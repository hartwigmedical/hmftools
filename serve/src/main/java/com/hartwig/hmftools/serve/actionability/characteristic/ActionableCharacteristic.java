package com.hartwig.hmftools.serve.actionability.characteristic;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableCharacteristic implements ActionableEvent {

    @NotNull
    public abstract TumorCharacteristicAnnotation name();

    @NotNull
    public abstract String cutOff();
}