package com.hartwig.hmftools.patientreporter.variants;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.Hotspot;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableSomaticVariant implements AllelicDepth {

    @NotNull
    public abstract String gene();

    public abstract boolean isDrupActionable();

    @NotNull
    public abstract String hgvsCodingImpact();

    @NotNull
    public abstract String hgvsProteinImpact();

    @Override
    public abstract int totalReadCount();

    @Override
    public abstract int alleleReadCount();

    @NotNull
    public abstract Hotspot hotspot();

    @NotNull
    public abstract Clonality clonality();

    @Nullable
    public abstract DriverCategory driverCategory();

    @Nullable
    public abstract Double driverLikelihood();

    public abstract double adjustedVAF();

    public abstract double adjustedCopyNumber();

    public abstract double minorAllelePloidy();

    public abstract boolean biallelic();
}
