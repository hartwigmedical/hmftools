package com.hartwig.hmftools.common.purple.purity;

import com.hartwig.hmftools.common.purple.PurpleQC;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurityQCContext {

    @NotNull
    public abstract PurityContext purityContext();

    @NotNull
    public abstract PurpleQC qc();
}
