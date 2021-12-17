package com.hartwig.hmftools.common.genome.position;

import org.jetbrains.annotations.Nullable;

public interface GenomeInterval extends GenomePosition {

    @Nullable
    Integer startOffset();

    default int startPosition() {
        Integer startOffset = startOffset();
        return startOffset == null ? position() : position() + startOffset;
    }

    @Nullable
    Integer endOffset();

    default int endPosition() {
        Integer endOffset = endOffset();
        return endOffset == null ? position() : position() + endOffset;
    }
}
