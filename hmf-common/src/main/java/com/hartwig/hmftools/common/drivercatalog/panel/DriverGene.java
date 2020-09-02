package com.hartwig.hmftools.common.drivercatalog.panel;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverGene extends Comparable<DriverGene> {

    @NotNull
    String gene();

    boolean reportMissenseAndInframe();

    boolean reportNonsenseAndFrameshift();

    boolean reportSplice();

    boolean reportDeletion();

    boolean reportDisruption();

    boolean reportAmplification();

    boolean reportHotspot();

    @NotNull
    DriverCategory likelihoodType();

    @Override
    default int compareTo(@NotNull final DriverGene o) {
        return gene().compareTo(o.gene());
    }

    default boolean reportVariant() {
        return reportMissenseAndInframe() || reportNonsenseAndFrameshift() || reportSplice() || reportHotspot();
    }
}
