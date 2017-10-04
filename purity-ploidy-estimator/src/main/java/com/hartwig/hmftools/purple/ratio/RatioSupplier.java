package com.hartwig.hmftools.purple.ratio;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.ReadRatio;

import org.jetbrains.annotations.NotNull;

public interface RatioSupplier {

    @NotNull
    Multimap<String, ReadRatio> tumorRatios();

    @NotNull
    Multimap<String, ReadRatio> referenceRatios();

}
