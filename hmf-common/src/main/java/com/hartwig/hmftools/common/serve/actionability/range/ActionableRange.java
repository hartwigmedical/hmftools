package com.hartwig.hmftools.common.serve.actionability.range;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.datamodel.range.RangeAnnotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableRange implements ActionableEvent, RangeAnnotation {

    @NotNull
    public abstract RangeType rangeType();
}
