package com.hartwig.hmftools.common.driver.panel;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;

import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCategory;

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
    boolean reportHetDeletion();
    boolean reportLoh();
    double hetDeletionThreshold();

    boolean reportDisruption();

    boolean reportAmplification();
    double amplificationRatio();

    boolean reportSomaticHotspot();

    @NotNull
    DriverGeneGermlineReporting reportGermlineVariant();

    @NotNull
    DriverGeneGermlineReporting reportGermlineHotspot();

    @NotNull
    DriverCategory likelihoodType();

    @NotNull
    DriverGeneGermlineReporting reportGermlineDeletion();

    @NotNull
    DriverGeneGermlineReporting reportGermlineDisruption();

    List<String> additionalReportedTranscripts();

    boolean reportPGX();

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
        return reportGermlineVariant() != NONE || reportGermlineHotspot() != NONE || reportGermlineDeletion() != NONE
                || reportGermlineDisruption() != NONE;
    }
}
