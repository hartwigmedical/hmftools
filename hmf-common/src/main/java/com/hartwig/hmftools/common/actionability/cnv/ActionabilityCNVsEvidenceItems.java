package com.hartwig.hmftools.common.actionability.cnv;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ActionabilityCNVsEvidenceItems {

    @NotNull
    public abstract List<ActionabilityCNVs> onLabel();

    @NotNull
    public abstract List<ActionabilityCNVs> offLabel();
}
