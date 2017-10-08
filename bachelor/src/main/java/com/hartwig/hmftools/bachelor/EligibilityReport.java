package com.hartwig.hmftools.bachelor;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class EligibilityReport {
    @NotNull
    public abstract String tag();

    @NotNull
    public abstract String program();

    @NotNull
    public abstract List<VariantContext> variants();
}