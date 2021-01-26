package com.hartwig.hmftools.ckb.variant;

import java.util.List;

import com.hartwig.hmftools.ckb.common.ReferenceInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantGeneVariantDescription {

    @NotNull
    public abstract String description();

    @NotNull
    public abstract List<ReferenceInfo> reference();


}
