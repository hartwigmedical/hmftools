package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurityContext {

    public abstract String version();

    public abstract Gender gender();

    @NotNull
    public abstract RunMode runMode();

    public abstract boolean targeted();

    public abstract FittedPurity bestFit();

    public abstract FittedPurityMethod method();

    public abstract FittedPurityScore score();

    public abstract PurpleQC qc();

    public abstract double polyClonalProportion();

    public abstract boolean wholeGenomeDuplication();

    public abstract double microsatelliteIndelsPerMb();

    public abstract double tumorMutationalBurdenPerMb();

    public abstract int tumorMutationalLoad();

    public abstract int svTumorMutationalBurden();

    @NotNull
    public abstract MicrosatelliteStatus microsatelliteStatus();

    @NotNull
    public abstract TumorMutationalStatus tumorMutationalLoadStatus();

    @NotNull
    public abstract TumorMutationalStatus tumorMutationalBurdenStatus();


    public static boolean checkHasReliablePurity(final PurityContext purityContext)
    {
        return purityContext.method() != FittedPurityMethod.NO_TUMOR;
    }
}
