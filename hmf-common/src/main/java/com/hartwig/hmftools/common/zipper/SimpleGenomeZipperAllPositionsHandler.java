package com.hartwig.hmftools.common.zipper;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface SimpleGenomeZipperAllPositionsHandler<R extends GenomeRegion, P extends GenomePosition> {
    void handle(@Nullable R region, @NotNull P position);
}
