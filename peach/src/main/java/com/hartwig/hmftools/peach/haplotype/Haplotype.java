package com.hartwig.hmftools.peach.haplotype;

import org.jetbrains.annotations.NotNull;

public interface Haplotype
{
    @NotNull
    String getName();

    boolean isWildType();
}
