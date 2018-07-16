package com.hartwig.hmftools.common.variant.cosmic;

import com.hartwig.hmftools.common.variant.TranscriptAnnotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CosmicAnnotation implements TranscriptAnnotation {

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String hgvsCoding();

    @NotNull
    public abstract String hgvsProtein();

    public abstract int count();

}
