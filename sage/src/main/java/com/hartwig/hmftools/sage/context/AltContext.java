package com.hartwig.hmftools.sage.context;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

public interface AltContext extends VariantHotspot {

    int rawSupportAlt();

    int rawDepth();

    int rawBaseQualityAlt();

    @NotNull
    String sample();

    @NotNull
    ReadContextCounter primaryReadContext();
}