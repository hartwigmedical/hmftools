package com.hartwig.hmftools.svanalysis.svgraph;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;

import java.util.Collection;
import java.util.List;

/**
 * Adjacency between DNA segments
 */
public interface BgAdjacency {
    BgSegment fromSegment();
    int fromOrientation();
    BgSegment toSegment();
    int toOrientation();
    BgEdge edge();
    Collection<String> linkedBy();
    default boolean isReference() {
        return edge().sv() == null;
    }
}
