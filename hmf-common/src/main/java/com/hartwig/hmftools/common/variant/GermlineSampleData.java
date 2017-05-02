package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
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
