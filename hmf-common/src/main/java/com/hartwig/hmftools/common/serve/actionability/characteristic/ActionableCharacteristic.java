package com.hartwig.hmftools.common.serve.actionability.characteristic;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.datamodel.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.common.serve.datamodel.characteristic.TumorCharacteristicsComparator;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableCharacteristic implements ActionableEvent {

    @NotNull
    public abstract TumorCharacteristicAnnotation name();

    @Nullable
    public abstract TumorCharacteristicsComparator comparator();

    @Nullable
    public abstract Double cutoff();
}