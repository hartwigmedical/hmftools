package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineSampleData {

    @NotNull
    public abstract String genoType();

    public abstract int totalReadCount();

    public abstract int alleleReadCount();

    public abstract int combinedDepth();

    public double alleleFrequency() {
        return (double) alleleReadCount() / totalReadCount();
    }
}
