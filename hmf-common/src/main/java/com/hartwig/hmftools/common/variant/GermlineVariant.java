package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineVariant implements Variant {
    public abstract static class Builder implements VariantBuilder<GermlineVariant> {}

    @NotNull
    public abstract GermlineSampleData refData();

    @Nullable
    public abstract GermlineSampleData tumorData();
}

