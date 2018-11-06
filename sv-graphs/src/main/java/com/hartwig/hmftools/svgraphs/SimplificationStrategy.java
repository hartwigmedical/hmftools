package com.hartwig.hmftools.svgraphs;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

import java.util.List;

public interface SimplificationStrategy {
    boolean shouldSimplify(Simplification s);
    boolean couldBeDirectlyLinked(EnrichedStructuralVariant sv, BgSegment segment);
    boolean couldBeDirectlyLinked(EnrichedStructuralVariant sv1, EnrichedStructuralVariant sv2);
    boolean couldBeFoldBackLinked(EnrichedStructuralVariant sv1, List<BgSegment> segment);
    boolean couldBeFoldBackLinked(EnrichedStructuralVariant sv1, EnrichedStructuralVariant sv2);
}
