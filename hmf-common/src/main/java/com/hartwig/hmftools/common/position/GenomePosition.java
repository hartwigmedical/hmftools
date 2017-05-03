package com.hartwig.hmftools.common.position;

import org.jetbrains.annotations.NotNull;

public interface GenomePosition {

    @NotNull
    String chromosome();

    long position();
}
