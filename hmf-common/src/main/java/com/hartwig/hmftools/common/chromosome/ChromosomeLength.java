package com.hartwig.hmftools.common.chromosome;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ChromosomeLength {

    @NotNull
    String chromosome();

    long length();
}
