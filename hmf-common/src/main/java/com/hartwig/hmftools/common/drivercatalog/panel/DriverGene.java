package com.hartwig.hmftools.common.drivercatalog.panel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverGene extends Comparable<DriverGene> {

    @NotNull
    String gene();

    @Nullable
    String deletionBand();

    boolean reportMissense();

    boolean reportTruncation();

    boolean reportSplice();

    boolean reportDisruption();

    boolean reportAmplification();

    boolean favorMultiHitAndBiallelic();

    @Override
    default int compareTo(@NotNull final DriverGene o) {
        return gene().compareTo(o.gene());
    }
}
