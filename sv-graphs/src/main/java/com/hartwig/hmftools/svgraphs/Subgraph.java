package com.hartwig.hmftools.svgraphs;

import java.util.Deque;
import java.util.List;
import java.util.Set;

public class Subgraph {
    public Set<BgAdjacency> externalConnections;
    public Set<BgAdjacency> internalConnections;
    public Set<BgSegment> internalSegments;
    public void calculatePathSets() {
    }
    public Deque<BgAdjacency> path;
    public List<BgAdjacency> paths;
}
