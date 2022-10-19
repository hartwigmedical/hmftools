package com.hartwig.hmftools.serve.sources.iclusion;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableTrial implements ActionableEvent {

}
