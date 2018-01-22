package com.hartwig.hmftools.common.cosmic.fusions;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class COSMICPromiscuousGene {
    @NotNull
    public abstract String geneName();

    @Nullable
    public abstract String transcript();
}
