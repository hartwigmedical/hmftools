package com.hartwig.hmftools.ckb.datamodelinterpretation.variant;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CategoryVariantPath {

    @NotNull
    public abstract String variantPath();

    @NotNull
    public abstract List<VariantInfo> variantInfo();
}
