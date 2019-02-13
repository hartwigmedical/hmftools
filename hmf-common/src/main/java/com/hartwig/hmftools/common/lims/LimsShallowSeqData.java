package com.hartwig.hmftools.common.lims;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class LimsShallowSeqData {

    @NotNull
    abstract String sampleId();

    abstract String purityShallowSeq();
}
