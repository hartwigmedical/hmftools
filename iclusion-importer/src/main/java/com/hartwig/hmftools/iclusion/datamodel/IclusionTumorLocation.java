package com.hartwig.hmftools.iclusion.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionTumorLocation {

    @NotNull
    public abstract String primaryTumorLocation();

    @NotNull
    public abstract List<String> doids();
}
