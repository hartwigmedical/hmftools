package com.hartwig.hmftools.serve.sources.actin;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActinTrial implements ActionableEvent {


}


