package com.hartwig.hmftools.common.position;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface GenomeInterval extends GenomePosition {

    default long startPosition() {
        return startOffset() == null ? position() : position() + startOffset();
    }

    default long endPosition() {
        return endOffset() == null ? position() : position() + endOffset();
    }

    @Nullable
    Integer startOffset();

    @Nullable
    Integer endOffset();

    default long intervalWidth() {
        return endPosition() - startPosition() + 1;
    }

    @NotNull
    default String chromosomeInterval() {
        return chromosome() + ":" + startPosition() + "-" + endPosition();
    }
}
