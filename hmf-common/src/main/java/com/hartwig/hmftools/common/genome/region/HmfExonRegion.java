package com.hartwig.hmftools.common.genome.region;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfExonRegion implements GenomeRegion {

    @NotNull
    public abstract int exonRank();
}
