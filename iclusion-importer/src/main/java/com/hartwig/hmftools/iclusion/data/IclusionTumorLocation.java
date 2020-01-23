package com.hartwig.hmftools.iclusion.data;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionTumorLocation {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String parentId();

    @NotNull
    public abstract String doid();

    @NotNull
    public abstract String doid2();

    @NotNull
    public abstract String indicationName();

    @NotNull
    public abstract String indicationNameFull();

    @NotNull
    public abstract List<String> nodeIds();
}
