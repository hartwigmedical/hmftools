package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.position.GenomePosition;

public interface BgEdge {
    Double ploidy();
    GenomePosition first();
    GenomePosition second();
}
