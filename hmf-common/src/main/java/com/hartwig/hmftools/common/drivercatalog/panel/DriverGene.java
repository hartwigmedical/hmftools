package com.hartwig.hmftools.common.drivercatalog.panel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverGene extends Comparable<DriverGene> {

    @NotNull
    String gene();

    @NotNull
    String deletionBand();

    boolean reportMissenseAndInframe();

    boolean reportNonsenseAndFrameshift();

    boolean reportSplice();

    boolean reportDeletionAndDisruption();

    boolean reportAmplification();

    boolean reportPromoterHotspots();

    @NotNull
    DriverLikelihoodType likelihoodType();

    @Override
    default int compareTo(@NotNull final DriverGene o) {
        return gene().compareTo(o.gene());
    }
}
