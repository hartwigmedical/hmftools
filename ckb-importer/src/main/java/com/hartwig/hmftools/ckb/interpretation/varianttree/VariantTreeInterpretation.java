package com.hartwig.hmftools.ckb.interpretation.varianttree;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantTreeInterpretation {

    @NotNull
    public abstract VariantInterpretation variantInterpretation();

    @NotNull
    public abstract GeneInterpretation geneInterpretation();
}
