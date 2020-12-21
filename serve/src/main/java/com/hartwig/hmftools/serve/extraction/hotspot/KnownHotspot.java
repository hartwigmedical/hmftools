package com.hartwig.hmftools.serve.extraction.hotspot;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.KnownEvent;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownHotspot implements VariantHotspot, KnownEvent {

    @NotNull
    public abstract String gene();

    @Nullable
    public abstract String transcript();

    @NotNull
    public abstract String proteinAnnotation();

}
