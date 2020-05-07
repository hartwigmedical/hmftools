package com.hartwig.hmftools.patientreporter.purple;

import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleSignatures {

    public abstract double microsatelliteIndelsPerMb();

    public abstract double tumorMutationalBurdenPerMb();

    public abstract int tumorMutationalLoad();

    @NotNull
    public abstract MicrosatelliteStatus microsatelliteStatus();

    @NotNull
    public abstract TumorMutationalStatus tumorMutationalLoadStatus();
}
