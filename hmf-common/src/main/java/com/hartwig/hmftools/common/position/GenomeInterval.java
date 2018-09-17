package com.hartwig.hmftools.common.position;

import org.jetbrains.annotations.Nullable;

public interface GenomeInterval extends GenomePosition {

    @Nullable
    Integer startOffset();

    default long startPosition() {
        Integer startOffset = startOffset();
        return startOffset == null ? position() : position() + startOffset;
    }

    @Nullable
    Integer endOffset();

    default long endPosition() {
        Integer endOffset = endOffset();
        return endOffset == null ? position() : position() + endOffset;
    }
}
