package com.hartwig.hmftools.common.sv;

import com.hartwig.hmftools.common.genome.position.GenomeInterval;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface StructuralVariantLeg extends GenomeInterval {

    byte orientation();

    @NotNull
    String homology();

    @Nullable
    Double alleleFrequency();

    @Nullable
    Integer inexactHomologyOffsetStart();

    @Nullable
    Integer inexactHomologyOffsetEnd();

    @Nullable
    Integer tumorVariantFragmentCount();

    @Nullable
    Integer tumorReferenceFragmentCount();

    @Nullable
    Integer normalVariantFragmentCount();

    @Nullable
    Integer normalReferenceFragmentCount();

    int anchoringSupportDistance();

    default int cnaPosition() {
        return orientation() ==  -1 ? position() : position() + 1;
    }
}