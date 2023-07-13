package com.hartwig.hmftools.datamodel.purple;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PurpleVariantTranscriptImpact {

    @NotNull
    String transcript();

    @NotNull
    String hgvsCodingImpact();

    @NotNull
    String hgvsProteinImpact();

    List<PurpleVariantEffect> effects();
}
