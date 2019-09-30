package com.hartwig.hmftools.sage.evidence;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SampleEvidence extends VariantHotspot {

    @NotNull
    String sample();

    int readDepth();

    int refSupport();

    int refQuality();

    int altSupport();

    int altQuality();

    int subprimeReadDepth();

    int readContextFull();

    int readContextPartial();

    int readContextRealigned();

    default double vaf() {
        return readDepth() == 0 ? 0 : (double) altSupport() / readDepth();
    }
}
