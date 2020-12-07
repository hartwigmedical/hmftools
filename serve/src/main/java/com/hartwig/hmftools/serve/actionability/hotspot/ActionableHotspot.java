package com.hartwig.hmftools.serve.actionability.hotspot;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionableHotspot implements ActionableEvent, VariantHotspot {

}


