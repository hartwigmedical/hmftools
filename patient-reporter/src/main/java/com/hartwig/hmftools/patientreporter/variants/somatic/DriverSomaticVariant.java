package com.hartwig.hmftools.patientreporter.variants.somatic;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DriverSomaticVariant {

    @NotNull
    public abstract SomaticVariant variant();

    public abstract double driverLikelihood();
}
