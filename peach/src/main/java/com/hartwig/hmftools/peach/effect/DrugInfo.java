package com.hartwig.hmftools.peach.effect;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrugInfo
{
    public abstract String drugName();
    public abstract String geneName();
    public abstract String prescriptionInfoUrl();
}
