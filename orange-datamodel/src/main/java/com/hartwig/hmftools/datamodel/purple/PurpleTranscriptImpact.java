package com.hartwig.hmftools.datamodel.purple;

import com.hartwig.hmftools.datamodel.variant.CodingEffect;
import com.hartwig.hmftools.datamodel.variant.impact.VariantEffect;
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
    public abstract Set<VariantEffect> effects();

    @NotNull
    public abstract CodingEffect codingEffect();

}
