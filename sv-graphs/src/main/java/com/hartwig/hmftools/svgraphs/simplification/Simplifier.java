package com.hartwig.hmftools.svgraphs.simplification;

import com.hartwig.hmftools.svgraphs.BreakpointGraph;

import java.util.ArrayList;
import java.util.List;

public abstract class Simplifier {
    protected final BreakpointGraph graph;
    protected final SimplificationStrategy strategy;
    public List<Simplification> findPotentialSimplifications() {
        return new ArrayList<>();
    }
    public void simplify(Simplification simplification) {
        throw new RuntimeException("NYI");
    }
    public Simplifier(BreakpointGraph graph, SimplificationStrategy strategy) {
        this.graph = graph;
        this.strategy = strategy;
    }
    public List<Simplification> simplify() {
        List<Simplification> simplified = new ArrayList<>();
        List<Simplification> list = findPotentialSimplifications();
        while (list.size() > 0) {
            for (Simplification s : list) {
                // might be invalidated by a previous simplification
                if (isValid(s)) {
                    simplify(s);
                    graph.sanityCheck();
                    simplified.add(s);
                }
            }
            graph.mergeReferenceSegments(false);
            list = findPotentialSimplifications();
        }
        return simplified;
    }

    /**
     * Determines whether the given simplication can be made based on the current state of the graph
     *
     * @param simplification
     * @return true if the simplification is valid, false if the given simplification is no longer possible
     */
    public boolean isValid(Simplification simplification) {
        return simplification.adjacencies()
                .stream()
                .allMatch(adj -> graph.getAllSegments().contains(adj.fromSegment()) && graph.getAllSegments().contains(adj.toSegment()) && graph.getOutgoing(
                        adj.fromSegment(),
                        adj.fromOrientation()).contains(adj));
    }
}
