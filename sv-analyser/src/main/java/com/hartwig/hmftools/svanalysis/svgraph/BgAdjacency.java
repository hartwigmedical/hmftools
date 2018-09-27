package com.hartwig.hmftools.svanalysis.svgraph;

/**
 * Adjacency between DNA segments
 */
public interface BgAdjacency {
    BgSegment fromSegment();
    int fromOrientation();
    BgSegment toSegment();
    int toOrientation();
    BgEdge edge();
}
