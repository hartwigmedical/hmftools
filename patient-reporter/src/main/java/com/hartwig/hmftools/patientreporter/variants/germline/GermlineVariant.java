package com.hartwig.hmftools.patientreporter.variants.germline;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.AllelicDepth;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineVariant implements GenomePosition, AllelicDepth {

    public abstract boolean passFilter();

    @NotNull
    public abstract String gene();

    @Override
    @NotNull
    public abstract String chromosome();

    @Override
    public abstract long position();

    @NotNull
    public abstract String hgvsCodingImpact();

    @NotNull
    public abstract String hgvsProteinImpact();

    @Override
    public abstract int totalReadCount();

    @Override
    public abstract int alleleReadCount();

    public abstract double adjustedVAF();

    public abstract double adjustedCopyNumber();

    public abstract boolean biallelic();

}
