package com.hartwig.hmftools.serve.sources.docm;

import com.hartwig.hmftools.serve.hotspot.HotspotSourceEntry;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DocmEntry implements HotspotSourceEntry {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String proteinAnnotation();

}
