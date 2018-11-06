package com.hartwig.hmftools.svgraphs;

import java.util.Collection;

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
