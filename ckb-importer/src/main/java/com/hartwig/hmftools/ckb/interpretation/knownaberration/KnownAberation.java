package com.hartwig.hmftools.ckb.interpretation.knownaberration;

import com.hartwig.hmftools.ckb.interpretation.common.variantinterpretation.VariantInterpretation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class KnownAberation {

    @NotNull
    public abstract VariantInterpretation knownAberations();
}
