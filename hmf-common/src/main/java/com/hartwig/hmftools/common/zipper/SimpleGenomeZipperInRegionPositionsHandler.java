package com.hartwig.hmftools.common.zipper;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@Deprecated
public interface SimpleGenomeZipperInRegionPositionsHandler<R extends GenomeRegion, P extends GenomePosition> {
    void handle(@NotNull R region, @NotNull P position);
}
