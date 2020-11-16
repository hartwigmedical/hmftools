package com.hartwig.hmftools.serve.hotspot;

import java.util.Set;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownHotspot implements VariantHotspot {

    @NotNull
    public abstract Set<String> sources();

    @NotNull
    public abstract String gene();

    @Nullable
    public abstract String transcript();

    @NotNull
    public abstract String proteinAnnotation();

}
