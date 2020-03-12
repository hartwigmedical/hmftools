package com.hartwig.hmftools.sage.context;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

public interface AltContextI extends VariantHotspot {

    int rawSupportRef();

    int rawSupportAlt();

    int rawDepth();

    double rawVaf();

    int rawBaseQualityRef();

    int rawBaseQualityAlt();

    @NotNull
    String sample();

}
