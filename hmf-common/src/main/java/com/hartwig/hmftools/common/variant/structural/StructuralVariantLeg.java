package com.hartwig.hmftools.common.variant.structural;

import com.hartwig.hmftools.common.position.GenomeInterval;

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
    Integer tumourVariantFragmentCount();

    @Nullable
    Integer tumourReferenceFragmentCount();

    @Nullable
    Integer normalVariantFragmentCount();

    @Nullable
    Integer normalReferenceFragmentCount();

    default long cnaPosition() {
        return orientation() ==  -1 ? position() : position() + 1;
    }
}