package com.hartwig.hmftools.orange.report.datamodel;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BreakendEntry {


    @NotNull
    public abstract String location();

    @NotNull
    public abstract String gene();

    public abstract boolean canonical();

    public abstract int exonUp();

    @NotNull
    public abstract StructuralVariantType type();

    @NotNull
    public abstract String range();

    public abstract int clusterId();

    public abstract double junctionCopyNumber();

    public abstract double undisruptedCopyNumber();

}
