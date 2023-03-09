package com.hartwig.hmftools.datamodel.purple;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Set;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleTranscriptImpact {

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String hgvsCodingImpact();

    @NotNull
    public abstract String hgvsProteinImpact();

    @Nullable
    public abstract Integer affectedCodon();

    @Nullable
    public abstract Integer affectedExon();

    public abstract boolean spliceRegion();

    @NotNull
    public abstract Set<PurpleVariantEffect> effects();

    @NotNull
    public abstract PurpleCodingEffect codingEffect();
}
