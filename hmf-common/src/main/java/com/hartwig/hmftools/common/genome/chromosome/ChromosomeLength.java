package com.hartwig.hmftools.common.genome.chromosome;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ChromosomeLength {

    @NotNull
    String chromosome();

    int length();
}
