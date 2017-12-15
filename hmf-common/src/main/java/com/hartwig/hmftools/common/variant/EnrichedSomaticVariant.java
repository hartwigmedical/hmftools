package com.hartwig.hmftools.common.variant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedSomaticVariant implements PurityAdjustedSomaticVariant {

    public abstract static class Builder implements PurityAdjustedSomaticVariantBuilder {
    }

    public abstract String trinucleotideContext();

    public abstract boolean highConfidenceRegion();

    public abstract int totalReadCount();

    public abstract int alleleReadCount();

    public abstract String microhomology();

    public abstract String refGenomeContext();

    public abstract String repeatSequence();

    public abstract String gene();

    public abstract String cosmicId();

    public abstract String dbsnpId();

    public abstract String effect();

    public abstract int repeatCount();

    public abstract double mappability();

    @Override
    public abstract double adjustedCopyNumber();

    @Override
    public abstract double adjustedVAF();
}
