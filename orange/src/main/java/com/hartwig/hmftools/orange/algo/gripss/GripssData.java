package com.hartwig.hmftools.orange.algo.gripss;

import java.util.List;

import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GripssData
{
    @NotNull
    public abstract List<StructuralVariant> allSomaticStructuralVariants();
}
