package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

/**
 * Adjacency between DNA segments
 */
public interface BgAdjacency {
    BgSegment fromSegment();
    int fromOrientation();
    BgSegment toSegment();
    int toOrientation();
    BgEdge edge();
    default boolean isReference() {
        return edge().sv() == null;
    }
}
