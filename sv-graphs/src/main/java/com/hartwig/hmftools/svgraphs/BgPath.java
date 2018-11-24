package com.hartwig.hmftools.svanalysis.svgraph;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Deque;
import java.util.List;

public class BgPath {
    private final BreakpointGraph bg;
    private final Deque<PathNode> path;

    public BgPath(BreakpointGraph bg, BgAdjacency seed) {
        this.bg = bg;
        this.path = new ArrayDeque<PathNode>();
        this.path.addLast(new PathNode(seed, seed.edge().ploidy()));
    }
    public static List<Collection<BgAdjacency>> potentialPaths(BreakpointGraph bg, BgAdjacency seed) {
        Deque<PathNode> path = new ArrayDeque<>();
        path.addLast(new PathNode(seed, seed.edge().ploidy()));
        // Path segments:
        // - events all have same CN; can traverse to end or cycle
        // - event matched with set of lower CN events
        // - double minute & fold-back inversions only possible way to increase CN without
        // adding teleomeres
    }
    public Collection<BgAdjacency> dfsRecurse() {
        for (BgAdjacency next : bg.nextBreakpointCandidates(path.getLast().edge)) {

        }
    }
    public static class PathNode {
        public PathNode(BgAdjacency edge, double copyNumber) {
            this.edge = edge;
            this.copyNumber = copyNumber;
        }
        BgAdjacency edge;
        double copyNumber;
        public boolean isFoldbackInversion() {
            return !edge.isReference() && edge.edge().sv().isFoldBackInversion();
        }
    }
}
