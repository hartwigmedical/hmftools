package com.hartwig.hmftools.common.actionability.somaticvariant;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ActionabilityRangeEvidenceItem {

    @NotNull
    public abstract List<ActionabilityRange> onLabel();

    @NotNull
    public abstract List<ActionabilityRange> offLabel();
}
