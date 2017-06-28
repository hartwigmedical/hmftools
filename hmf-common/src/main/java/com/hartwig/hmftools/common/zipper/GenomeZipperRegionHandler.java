package com.hartwig.hmftools.common.zipper;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public interface GenomeZipperRegionHandler<R extends GenomeRegion> {

    void chromosome(@NotNull String chromosome);

    void enter(@NotNull R region);

    void exit(@NotNull R region);
}