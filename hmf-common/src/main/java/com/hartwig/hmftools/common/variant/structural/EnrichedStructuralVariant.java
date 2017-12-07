package com.hartwig.hmftools.common.variant.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedStructuralVariant implements StructuralVariant {

    @Nullable
    public abstract Double ploidy();

    @Nullable
    public abstract Double adjustedStartAF();

    @Nullable
    public abstract Double adjustedStartCopyNumber();

    @Nullable
    public abstract Double adjustedStartCopyNumberChange();

    @Nullable
    public abstract Double adjustedEndAF();

    @Nullable
    public abstract Double adjustedEndCopyNumber();

    @Nullable
    public abstract Double adjustedEndCopyNumberChange();
}
