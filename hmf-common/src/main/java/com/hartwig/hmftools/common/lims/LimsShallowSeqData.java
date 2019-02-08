package com.hartwig.hmftools.common.lims;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class LimsShallowSeqData {
    // only known if shallow seq is done

    @Nullable
    abstract String sampleId();

    @Nullable
    abstract String purityShallowSeq();
}
