package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedSomaticVariant implements Variant {

    public abstract boolean highConfidenceRegion();

    public abstract int totalReadCount();

    public abstract int alleleReadCount();

    public double alleleFrequency() {
        return (double) alleleReadCount() / totalReadCount();
    }

    @Value.Default
    public double adjustedCopyNumber() {
        return 0;
    }

    @Value.Default
    public double adjustedVAF() {
        return 0;
    }
}
