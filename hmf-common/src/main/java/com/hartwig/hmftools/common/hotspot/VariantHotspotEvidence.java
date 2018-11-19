package com.hartwig.hmftools.common.hotspot;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VariantHotspotEvidence extends VariantHotspot {

    int readDepth();

    int refSupport();

    int altSupport();

    int altQuality();

    int indelSupport();
}
