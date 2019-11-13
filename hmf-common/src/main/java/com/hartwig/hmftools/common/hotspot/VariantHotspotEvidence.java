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

    int refQuality();

    int altSupport();

    int altQuality();

    int indelSupport();

    int altMapQuality();

    int altMinQuality();

    int altDistanceFromRecordStart();

    int altMinDistanceFromAlignment();

    int subprimeReadDepth();

    String readContext();

    int readContextCount();

    int readContextCountOther();

    default double vaf() {
        return readDepth() == 0 ? 0 : (double) altSupport() / readDepth();
    }
}
