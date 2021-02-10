package com.hartwig.hmftools.ckb.interpretation.common.variantinterpretation;

import com.hartwig.hmftools.ckb.datamodelinterpretation.gene.Gene;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantInterpretation {

    @Nullable
    public abstract Variant variant();

    @Nullable
    public abstract Gene gene();

}
