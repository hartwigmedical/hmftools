package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleCharacteristics {

    public abstract boolean wholeGenomeDuplication();

    public abstract double microsatelliteIndelsPerMb();

    @NotNull
    public abstract MicrosatelliteStatus microsatelliteStatus();

    public abstract double tumorMutationalBurdenPerMb();

    public abstract int tumorMutationalLoad();

    @NotNull
    public abstract TumorMutationalStatus tumorMutationalLoadStatus();

    public abstract int svTumorMutationalBurden();
}
