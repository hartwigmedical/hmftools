package com.hartwig.hmftools.common.position;

import org.jetbrains.annotations.NotNull;

public interface GenomeInterval extends GenomePosition {

    long startPosition();

    long endPosition();

    default long intervalWidth() {
        return endPosition() - startPosition() + 1;
    }

    @NotNull
    default String chromosomeInterval() {
        return chromosome() + ":" + startPosition() + "-" + endPosition();
    }
}
