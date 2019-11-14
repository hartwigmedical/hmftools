package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableVariant implements Variant {

    @NotNull
    @Override
    public abstract String gene();

    @NotNull
    @Override
    public abstract String chromosome();

    @Override
    public abstract long position();

    @NotNull
    @Override
    public abstract String ref();

    @NotNull
    @Override
    public abstract String alt();

    @NotNull
    @Override
    public abstract CodingEffect canonicalCodingEffect();

    @NotNull
    public abstract String gDNA();

    @NotNull
    public abstract String hgvsCodingImpact();

    @NotNull
    public abstract String hgvsProteinImpact();

    @Override
    public abstract int totalReadCount();

    @Override
    public abstract int alleleReadCount();

    public abstract double totalPloidy();

    public abstract double allelePloidy();

    @NotNull
    public abstract Hotspot hotspot();

    public abstract double clonalLikelihood();

    @Nullable
    public abstract DriverCategory driverCategory();

    @Nullable
    public abstract Double driverLikelihood();

    public abstract boolean biallelic();

    public abstract boolean notifyClinicalGeneticist();
}
