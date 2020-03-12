package com.hartwig.hmftools.sage.context;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

public interface AltContext extends VariantHotspot {

    int rawSupportRef();

    int rawSupportAlt();

    int rawDepth();

    double rawVaf();

    int rawBaseQualityRef();

    int rawBaseQualityAlt();

    @NotNull
    String sample();

    @NotNull
    ReadContextCounter primaryReadContext();
}