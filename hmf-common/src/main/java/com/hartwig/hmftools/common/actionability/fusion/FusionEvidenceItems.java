package com.hartwig.hmftools.common.actionability.fusion;

import java.util.List;

import com.hartwig.hmftools.common.actionability.somaticvariant.EvidenceItem;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class FusionEvidenceItems {

    @NotNull
    public abstract List<ActionabilityFusionPairs> onLabel();

    @NotNull
    public abstract List<ActionabilityFusionPairs> offLabel();

}
