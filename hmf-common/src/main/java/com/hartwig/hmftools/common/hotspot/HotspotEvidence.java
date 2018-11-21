package com.hartwig.hmftools.common.hotspot;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface HotspotEvidence extends GenomePosition {

    @NotNull
    HotspotEvidenceType type();

    @NotNull
    String ref();

    @NotNull
    String alt();

    int qualityScore();

    int tumorRefCount();

    int tumorAltCount();

    int tumorReads();

    int normalRefCount();

    int normalAltCount();

    int normalReads();

    int normalIndelCount();

    default boolean isIndel() {
        return ref().length() != alt().length();
    }
}
