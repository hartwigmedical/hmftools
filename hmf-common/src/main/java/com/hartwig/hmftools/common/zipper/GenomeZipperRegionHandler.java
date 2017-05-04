package com.hartwig.hmftools.common.zipper;

import com.hartwig.hmftools.common.region.GenomeRegion;

public interface GenomeZipperRegionHandler<R extends GenomeRegion> {
    void enter(R region);

    void exit(R region);
}