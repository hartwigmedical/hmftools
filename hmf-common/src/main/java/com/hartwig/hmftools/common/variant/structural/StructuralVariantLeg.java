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

    default Integer inexactHomologyLength() {
        return inexactHomologyOffsetStart() != null && inexactHomologyOffsetEnd() != null ?
                inexactHomologyOffsetEnd() - inexactHomologyOffsetStart() : null;
    }

    @Nullable
    Integer inexactHomologyOffsetStart();

    @Nullable
    Integer inexactHomologyOffsetEnd();

    default Integer tumourFragmentCount() {
        return tumourVariantFragmentCount() != null && tumourReferenceFragmentCount() != null ?
                tumourVariantFragmentCount() + tumourReferenceFragmentCount() : null;
    }

    default Integer normalFragmentCount() {
        return normalVariantFragmentCount() != null && normalReferenceFragmentCount() != null ?
                normalVariantFragmentCount() + normalReferenceFragmentCount() : null;
    }

    @Nullable
    Integer tumourVariantFragmentCount();

    @Nullable
    Integer tumourReferenceFragmentCount();

    @Nullable
    Integer normalVariantFragmentCount();

    @Nullable
    Integer normalReferenceFragmentCount();
}