package com.hartwig.hmftools.serve.sources.vicc;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ActionableEvidence implements ActionableEvent {

}
