package com.hartwig.hmftools.patientreporter.germline;

import com.hartwig.hmftools.common.variant.AllelicDepth;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineVariant implements AllelicDepth {

    public abstract boolean passFilter();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String hgvsCodingImpact();

    @NotNull
    public abstract String hgvsProteinImpact();

    @Override
    public abstract int totalReadCount();

    @Override
    public abstract int alleleReadCount();

    @NotNull
    public abstract String germlineStatus();

    public abstract double adjustedVAF();

    public abstract double adjustedCopyNumber();

    public abstract double minorAllelePloidy();

    public abstract boolean biallelic();

}
