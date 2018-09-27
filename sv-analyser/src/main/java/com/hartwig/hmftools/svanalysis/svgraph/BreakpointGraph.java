package com.hartwig.hmftools.svanalysis.svgraph;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;
import java.util.stream.Collectors;

public class BreakpointGraph {
    private static final Logger LOGGER = LogManager.getLogger(BreakpointGraph.class);

    public BgSegment UnplacedSegment = BgSegment.createUnplacedSegment();
    Map<BgSegment, List<BgAdjacency>> startEdges;
    Map<BgSegment, List<BgAdjacency>> endEdges;

    public BreakpointGraph(List<PurpleCopyNumber> cnRecords, List<EnrichedStructuralVariant> svRecords) {
        Map<String, List<BgSegment>> segments = cnRecords.stream()
                .map(cn -> new BgSegment(cn))
                .sorted(ByGenomicPosition)
                .collect(Collectors.groupingBy(GenomePosition::chromosome, Collectors.toList()));
        for (List<BgSegment> chrSegs: segments.values()) {
            for (int i = 0; i < chrSegs.size(); i++) {
                BgSegment s = chrSegs.get(i);
                startEdges.put(s, new ArrayList<>());
                endEdges.put(s, new ArrayList<>());
            }
            for (int i = 0; i < chrSegs.size() - 1; i++) {
                BgSegment left = chrSegs.get(i);
                BgSegment right = chrSegs.get(i + 1);
                BgReferenceEdge edge = new BgReferenceEdge(left, right);
                addEdge(left, 1, right, -1, edge);
            }
        }
        for (EnrichedStructuralVariant sv : svRecords) {
            addSV(sv, segments);
        }
    }
    private void addSV(EnrichedStructuralVariant sv, Map<String, List<BgSegment>> segmentLookup) {
        BgSegment startSegment = containingSegment(sv.start(), segmentLookup);
        BgSegment endSegment  = containingSegment(sv.end(), segmentLookup);
        // sanity checks
        if (sv.start().orientation() == 1) {
            assert(startSegment.endPosition() == sv.start().position());
        } else {
            assert(startSegment.startPosition() == sv.start().position());
        }
        if (sv.end() != null) {
            if (sv.end().orientation() == 1) {
                assert(endSegment.endPosition() == sv.end().position());
            } else {
                assert(endSegment.startPosition() == sv.end().position());
            }
        }
        addEdge(startSegment, startSegment == UnplacedSegment ? -1 :  sv.start().orientation(),
                endSegment, endSegment == UnplacedSegment ? -1 : sv.end().orientation(),
                new BgSv(sv));
    }
    private BgSegment containingSegment(EnrichedStructuralVariantLeg leg, Map<String, List<BgSegment>> segmentLookup) {
        if (leg == null) {
            return UnplacedSegment;
        }
        List<BgSegment> segments = segmentLookup.get(leg.chromosome());
        if (segments == null) {
            return UnplacedSegment;
        }
        int position = Collections.binarySearch(segments, leg);
        if (position < 0) {
            position = 2 - position;
        }
        if (position >= 0 && position < segments.size()) {
            return segments.get(position);
        }
        return UnplacedSegment;
    }
    private void addEdge(BgSegment left, int leftOrientation, BgSegment right, int rightOrientation, BgEdge edge) {
        BgAdjacency leftAdj = ImmutableBgAdjacencyImpl.builder()
                .fromSegment(left)
                .fromOrientation(leftOrientation)
                .toSegment(right)
                .toOrientation(rightOrientation)
                .edge(edge)
                .build();
        BgAdjacency rightAdj = ImmutableBgAdjacencyImpl.builder()
                .fromSegment(right)
                .fromOrientation(rightOrientation)
                .toSegment(left)
                .toOrientation(leftOrientation)
                .edge(edge)
                .build();
        addAdjacency(leftAdj);
        addAdjacency(rightAdj);
    }
    private void addAdjacency(BgAdjacency adj) {
        switch (adj.fromOrientation()) {
            case 1:
                endEdges.get(adj.fromSegment()).add(adj);
                break;
            case -1:
                startEdges.get(adj.fromSegment()).add(adj);
                break;
            default:
                throw new IllegalArgumentException("Invalid orientation");
        }
    }
    public List<StructuralVariant> simplify() {
        List<StructuralVariant> simplified = new ArrayList<>();
        StructuralVariant var = simplifySimpleDuplications();
        while (var != null) {
            simplified.add(var);
            mergeReferenceSegments();
            var = simplifySimpleDuplications();
        }
        return simplified;
    }
    public int mergeReferenceSegments() {
        for (BgSegment segment : endEdges.keySet()) {
            BgSegment nextSegment = nextReferenceSegment(segment);
            if (nextSegment != null && svCount(getOutgoing(segment, 1)) == 0 && svCount(getOutgoing(segment, -1)) == 0) {
                merge(segment, nextSegment);
                return 1 + mergeReferenceSegments();
            }
        }
        return 0;
    }
    public StructuralVariant simplifySimpleDuplications() {
        for (BgSegment segment : startEdges.keySet()) {
            List<BgAdjacency> prev = getOutgoing(segment, -1);
            List<BgAdjacency> next = getOutgoing(segment, 1);
            // simple duplication check
            if (svCount(prev) == 1 || svCount(next) == 1) {
                BgSv prevSv = firstSv(prev);
                BgSv nextSv = firstSv(next);
                if (prevSv != null && prevSv == nextSv) {
                    LOGGER.debug("Simple tandem duplication {}-{}", prevSv.first(), nextSv.second());
                    // TODO: check CN consistency
                    startEdges.get(segment).remove(prev);
                    endEdges.get(segment).remove(prev);
                    return prevSv.sv();
                }
            }
        }
        return null;
    }
    public BgSegment prevReferenceSegment(BgSegment segment) {
        for (BgAdjacency adj : getOutgoing(segment, -1)) {
            if (adj.edge() instanceof BgReferenceEdge) {
                return adj.toSegment();
            }
        }
        return null;
    }
    public BgSegment nextReferenceSegment(BgSegment segment) {
        for (BgAdjacency adj : getOutgoing(segment, 1)) {
            if (adj.edge() instanceof BgReferenceEdge) {
                return adj.toSegment();
            }
        }
        return null;
    }
    public static BgSv firstSv(List<BgAdjacency> adjacencies) {
        for (BgAdjacency adj : adjacencies) {
            if (adj.edge() instanceof BgSv) {
                return (BgSv)adj;
            }
        }
        return null;
    }
    public static int svCount(List<BgAdjacency> adjacencies) {
        int count = 0;
        for (BgAdjacency adj : adjacencies) {
            if (adj.edge() instanceof BgSv) {
                count++;
            }
        }
        return count;
    }
    public List<BgAdjacency> getOutgoing(BgSegment segment, int orientation) {
        switch (orientation) {
            case 1:
                return endEdges.get(segment);
            case -1:
                return startEdges.get(segment);
            default:
                throw new IllegalArgumentException("Invalid orientation");
        }
    }
    public BgSegment merge(BgSegment left, BgSegment right) {
        List<BgAdjacency> leftStart = startEdges.get(left);
        List<BgAdjacency> rightStart = startEdges.get(right);
        List<BgAdjacency> leftEnd = endEdges.get(left);
        List<BgAdjacency> rightEnd = endEdges.get(right);

        if (leftEnd.size() != 1
                || rightStart.size() == 1
                || !(leftEnd.get(0) instanceof BgReferenceEdge)
                || !(rightStart.get(0) instanceof BgReferenceEdge)) {
            throw new IllegalArgumentException("Cannot merge DNA segments separated by a SV");
        }
        // Replace segment
        BgSegment merged = BgSegment.merge(left, right);
        startEdges.remove(left);
        startEdges.remove(right);
        endEdges.remove(left);
        endEdges.remove(right);
        startEdges.put(merged, leftStart);
        endEdges.put(merged, rightEnd);
        return merged;
    }
    private static final Ordering<GenomePosition> ByGenomicPosition = new Ordering<GenomePosition>() {
        public int compare(GenomePosition o1, GenomePosition o2) {
            return ComparisonChain.start()
                    .compare(o1.chromosome(), o2.chromosome())
                    .compare(o1.position(), o2.position())
                    .result();
        }
    };
}
