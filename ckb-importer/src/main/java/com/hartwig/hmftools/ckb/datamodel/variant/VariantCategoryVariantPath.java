package com.hartwig.hmftools.ckb.datamodel.variant;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.VariantInfo;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantCategoryVariantPath {

    @NotNull
    public abstract String variantPath();

    @NotNull
    public abstract List<VariantInfo> variant();


}
