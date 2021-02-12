package com.hartwig.hmftools.ckb.datamodelinterpretation.variant;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodelinterpretation.common.ReferenceExtend;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantDescription {

    @Nullable
    public abstract String description();

    @NotNull
    public abstract List<ReferenceExtend> references();
}
