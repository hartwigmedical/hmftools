package com.hartwig.hmftools.common.drivercatalog.panel;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverGene extends Comparable<DriverGene>
{
    @NotNull
    String gene();

    boolean reportMissenseAndInframe();

    boolean reportNonsenseAndFrameshift();

    boolean reportSplice();

    boolean reportDeletion();

    boolean reportDisruption();

    boolean reportAmplification();

    boolean reportSomaticHotspot();

    @NotNull
    DriverGeneGermlineReporting reportGermlineVariant();

    @NotNull
    DriverGeneGermlineReporting reportGermlineHotspot();

    @NotNull
    DriverCategory likelihoodType();

    boolean reportGermlineDisruption();

    @Nullable
    String alternativeTranscripts();

    @Override
    default int compareTo(@NotNull final DriverGene o)
    {
        return gene().compareTo(o.gene());
    }

    default boolean reportSomatic()
    {
        return reportMissenseAndInframe() || reportNonsenseAndFrameshift() || reportSplice() || reportSomaticHotspot();
    }

    default boolean reportGermline()
    {
        return reportGermlineVariant() != DriverGeneGermlineReporting.NONE || reportGermlineHotspot() != DriverGeneGermlineReporting.NONE;
    }
}
