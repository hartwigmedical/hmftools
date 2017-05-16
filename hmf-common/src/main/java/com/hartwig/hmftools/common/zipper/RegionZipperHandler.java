package com.hartwig.hmftools.common.zipper;

import com.hartwig.hmftools.common.region.GenomeRegion;

public interface RegionZipperHandler<S extends GenomeRegion, T extends GenomeRegion> {
    void left(S region);
    void right(T region);
}