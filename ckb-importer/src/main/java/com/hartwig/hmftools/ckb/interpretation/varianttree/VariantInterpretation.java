package com.hartwig.hmftools.ckb.interpretation.varianttree;

import com.hartwig.hmftools.ckb.datamodel.gene.Gene;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantInterpretation {

    @NotNull
    public abstract Variant variant();

    @Nullable
    public abstract Gene gene();

}
