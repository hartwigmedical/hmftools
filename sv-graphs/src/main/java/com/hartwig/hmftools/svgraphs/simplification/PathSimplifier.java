package com.hartwig.hmftools.svgraphs.simplification;

import com.hartwig.hmftools.svgraphs.BgAdjacency;
import com.hartwig.hmftools.svgraphs.BreakpointGraph;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public class PathSimplifier extends Simplifier {
    private static final Logger LOGGER = LogManager.getLogger(PathSimplifier.class);

    public PathSimplifier(BreakpointGraph graph, SimplificationStrategy strategy) {
        super(graph, strategy);
    }
    public Simplification findPathTo(BgAdjacency from, BgAdjacency to) {
        throw new NotImplementedException();
    }
    /*
    // TODO: outstanding simplification
    // - Simple INS/translocation
    // - Simple inversion
    // - Fold-back reduction to next breakpoints (do we even need to special case this?)
    // - match highest CN breakend to partner
    public void findStronglyConnectedComponentHierarchy() {
    }
    public List<Set<BgPath>> candidateSubgraphResolutions() {
        // Do fold-back walks from highest ploidy event
        foldback_traversal();

        // Takes a subgraph and enumerates all possible resolution paths
        // Walk along, or down; consuming breakpoints as we go.
        // score the traversal set based on consistency.
        // Traverse greedily, step back

        // TODO: how do we identify which subgraph?
        //  - try clustering algorithms
        // TODO: initial implementation: brute force through paths
        Set<BgAdjacency> subgraphExternalLinks;
        // Options:
        // 1) start at an external link and walk in?
        // 2) start at the max CN event and expand in both directions
    }
    //endregion
    //region Traversals
    public void test_traverse(BgAdjacency adj) {
        Deque<BgAdjacency> path = new ArrayDeque<>();
        path.add(adj);
        traverse_recursive(path, adj.edge().ploidy());
    }
    private void traverse_recursive(Deque<BgAdjacency> path, double nominalCn) {
        BgAdjacency adj = getPartner(path.getLast());
        Collection<BgAdjacency> nextBp = Stream.concat(
                nextBreakpointCandidates(adj).stream(),
                nextFoldbackDoublingCandidates(adj).stream()).collect(Collectors.toList());
        for (BgAdjacency nextAdj : nextBp) {
            if (nextAdj == null) { // or outside of subgraph bounds
                traverse_completed(path);
            }
            // Adjust CN of intervening segments
            List<BgSegment> seg = getSegmentsBetween(adj, nextAdj);
            seg.stream().forEach(s -> s.pushCopyNumberChange(-nominalCn));
            path.addLast(adj); // or is this the partner?
            if (adj.toSegment() == getUnplacedSegment()) {
                traverse_completed(path);
            } else {
                traverse_recursive(path, nextAdj.edge().ploidy());
            }
            path.removeLast();
            seg.stream().forEach(s -> s.popCopyNumberChange(nominalCn));
        }
    }
    private void traverse_completed(Deque<BgAdjacency> path) {
    }
    */
}
