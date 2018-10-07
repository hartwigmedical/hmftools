package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

public interface BgEdge {
    Double ploidy();
    GenomePosition first();
    GenomePosition second();
    EnrichedStructuralVariant sv();
}
