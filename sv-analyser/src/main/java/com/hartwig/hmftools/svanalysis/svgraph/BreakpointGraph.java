package com.hartwig.hmftools.svanalysis.svgraph;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterables;
import com.google.common.collect.Ordering;
import com.google.common.collect.Streams;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.util.*;
import java.util.stream.Collectors;

public class BreakpointGraph {
    private static final Logger LOGGER = LogManager.getLogger(BreakpointGraph.class);

    private final BgSegment unplacedSegment = BgSegment.createUnplacedSegment();
    private final Map<BgSegment, List<BgAdjacency>> startEdges = new HashMap<>();
    private final Map<BgSegment, List<BgAdjacency>> endEdges = new HashMap<>();

    private BreakpointGraph() {
        startEdges.put(unplacedSegment, new ArrayList<>());
        endEdges.put(unplacedSegment, new ArrayList<>());
    }
    public BreakpointGraph(List<PurpleCopyNumber> cnRecords, List<EnrichedStructuralVariant> svRecords) {
        this();
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
        sanityCheck();
    }
    public BgSegment getUnplacedSegment() {
        return unplacedSegment;
    }
    public Collection<BgSegment> getAllSegments() {
        return startEdges.keySet();
    }
    public Collection<EnrichedStructuralVariant> getAllStructuralVariants() {
        return getAllEdges().stream()
                .filter(e -> e instanceof BgSv)
                .map(e -> ((BgSv)e).sv())
                .collect(Collectors.toList());
    }
    public Collection<BgEdge> getAllEdges() {
        Set<BgEdge> edges = Streams.concat(startEdges.values().stream(), endEdges.values().stream())
                .flatMap(adjList -> adjList.stream())
                .map(adj -> adj.edge())
                .collect(Collectors.toSet());
        return edges;
    }
    public BgSegment getSegment(PurpleCopyNumber purpleCopyNumber) {
        return startEdges.keySet().stream()
                .filter(s -> s.cn() == purpleCopyNumber)
                .findFirst()
                .orElse(null);
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
        addEdge(startSegment, startSegment == unplacedSegment ? -1 :  sv.start().orientation(),
                endSegment, endSegment == unplacedSegment ? -1 : sv.end().orientation(),
                new BgSv(sv));
    }
    private BgSegment containingSegment(EnrichedStructuralVariantLeg leg, Map<String, List<BgSegment>> segmentLookup) {
        if (leg == null) {
            return unplacedSegment;
        }
        List<BgSegment> segments = segmentLookup.get(leg.chromosome());
        if (segments == null) {
            return unplacedSegment;
        }
        int position = Collections.binarySearch(segments, leg);
        if (position < 0) {
            position = -2 - position;
        }
        if (position >= 0 && position < segments.size()) {
            return segments.get(position);
        }
        return unplacedSegment;
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
    public StructuralVariant simplifySimpleDeletion() {
        for (BgSegment segment : startEdges.keySet()) {
            BgSegment prevSegment = prevReferenceSegment(segment);
            BgSegment nextSegment = nextReferenceSegment(segment);
            if (prevSegment != null && nextSegment != null &&
                    svCount(getOutgoing(segment, -1)) == 0 &&
                    svCount(getOutgoing(segment, 1)) == 0 &&
                    svCount(getOutgoing(prevSegment, 1)) == 1 &&
                    svCount(getOutgoing(nextSegment, -1)) == 1) {
                BgAdjacency leftSv = firstSv(getOutgoing(prevSegment, 1));
                BgAdjacency rightSv = firstSv(getOutgoing(nextSegment, -1));
                if (leftSv.edge() == rightSv.edge()) {
                    // we have a deletion
                    if (shouldCollapseDeletion(prevSegment, segment, nextSegment, leftSv.edge())) {
                        LOGGER.debug("Simplifying simple deletion of {}", segment);
                        removeAdjacency(leftSv, rightSv);
                        // TODO: adjust CN of merged segment by the ploidy of the deletion
                        // we don't actually want a straight merge
                        segment = merge(prevSegment, segment);
                        segment = merge(segment, nextSegment);
                        return leftSv.edge().sv();
                    }
                }
            }
        }
        return null;
    }

    private boolean shouldCollapseDeletion(BgSegment prevSegment, BgSegment segment, BgSegment nextSegment, BgEdge edge) {
        // TODO: check for CN consistency
        return true;
    }
    private boolean shouldCollapseDuplication(BgSegment prevSegment, BgSegment segment, BgSegment nextSegment, BgEdge edge) {
        // TODO: check for CN consistency
        return true;
    }

    public EnrichedStructuralVariant simplifySimpleDuplications() {
        for (BgSegment segment : startEdges.keySet()) {
            List<BgAdjacency> prev = startEdges.get(segment);
            List<BgAdjacency> next = endEdges.get(segment);
            // simple duplication check
            if (svCount(prev) == 1 && svCount(next) == 1) {
                BgAdjacency prevSv = firstSv(prev);
                BgAdjacency nextSv = firstSv(next);
                if (prevSv.edge() == nextSv.edge()) {
                    if (shouldCollapseDuplication(prevReferenceSegment(segment), segment, nextReferenceSegment(segment), prevSv.edge())) {
                        LOGGER.debug("Simplifying simple tandem duplication of {}", segment);
                        removeAdjacency(prevSv, nextSv);
                        // TODO: adjust CN of segment by the ploidy of the duplication
                        return prevSv.edge().sv();
                    }
                }
            }
        }
        return null;
    }

    private void removeAdjacency(@NotNull BgAdjacency leftAdj, @NotNull BgAdjacency rightAdj) {
        if (leftAdj == null || rightAdj == null) {
            throw new NullPointerException();
        }
        if (leftAdj.edge() != rightAdj.edge()) {
            throw new IllegalArgumentException("Adjacencies are not paired.");
        }
        if (leftAdj == rightAdj) {
            throw new IllegalArgumentException("Adjacencies are the same object");
        }
        if (!getOutgoing(leftAdj.fromSegment(), leftAdj.fromOrientation()).remove(leftAdj)) {
            throw new IllegalStateException("Sanity check failure: removed non-existent adjacency");
        }
        if (!getOutgoing(rightAdj.fromSegment(), rightAdj.fromOrientation()).remove(rightAdj)) {
            throw new IllegalStateException("Sanity check failure: removed non-existent adjacency");
        }
        sanityCheck();
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
    public static BgAdjacency firstSv(List<BgAdjacency> adjacencies) {
        for (BgAdjacency adj : adjacencies) {
            if (!adj.isReference()) {
                return adj;
            }
        }
        return null;
    }
    public static int svCount(List<BgAdjacency> adjacencies) {
        int count = 0;
        for (BgAdjacency adj : adjacencies) {
            if (!adj.isReference()) {
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
        if (left == unplacedSegment || right == unplacedSegment) {
            throw new IllegalArgumentException("Cannot merge placeholder unplaced DNA segment");
        }
        List<BgAdjacency> leftStart = startEdges.get(left);
        List<BgAdjacency> rightStart = startEdges.get(right);
        List<BgAdjacency> leftEnd = endEdges.get(left);
        List<BgAdjacency> rightEnd = endEdges.get(right);

        if (svCount(leftEnd) > 0 || svCount(rightStart) > 0) {
            throw new IllegalArgumentException("Cannot merge DNA segments separated by a SV");
        }
        // Create new segment
        BgSegment merged = BgSegment.merge(left, right);
        // update adjacency from segments
        startEdges.put(merged, leftStart.stream()
            .map(adj -> ImmutableBgAdjacencyImpl.builder()
                .from(adj)
                .fromSegment(merged)
                .build())
            .collect(Collectors.toList()));
        endEdges.put(merged, rightEnd.stream()
            .map(adj -> ImmutableBgAdjacencyImpl.builder()
                    .from(adj)
                    .fromSegment(merged)
                    .build())
            .collect(Collectors.toList()));
        // update adjacency to segments
        for (BgAdjacency adj: Iterables.concat(leftStart, rightEnd)) {
            BgAdjacency partner = getPartner(adj);
            BgAdjacency newPartner = ImmutableBgAdjacencyImpl.builder()
                    .from(partner)
                    .toSegment(merged)
                    .build();
            getOutgoing(adj.toSegment(), adj.toOrientation()).set(getOutgoing(adj.toSegment(), adj.toOrientation()).indexOf(partner), newPartner);
        }
        // remove old segments from the graph
        startEdges.remove(left);
        startEdges.remove(right);
        endEdges.remove(left);
        endEdges.remove(right);
        LOGGER.debug("Merged CN: {}:{}-{} CN={} with {}:{}-{} CN={} to create {}:{}-{} CN={}",
                left.chromosome(), left.startPosition(), left.endPosition(), left.ploidy(),
                right.chromosome(), right.startPosition(), right.endPosition(), right.ploidy(),
                merged.chromosome(), merged.startPosition(), merged.endPosition(), merged.ploidy());
        sanityCheck();
        return merged;
    }
    private BgAdjacency getPartner(BgAdjacency adj) {
        for (BgAdjacency remoteAdj : getOutgoing(adj.toSegment(), adj.toOrientation())) {
            if (remoteAdj != adj && remoteAdj.edge() == adj.edge()) {
                return remoteAdj;
            }
        }
        throw new IllegalStateException("Partner adjacency missing from graph");
    }
    private void sanityCheck() {
        assert(startEdges.keySet().containsAll(endEdges.keySet()));
        assert(endEdges.keySet().containsAll(startEdges.keySet()));
        for (BgSegment segment : startEdges.keySet()) {
            for (BgAdjacency adj : startEdges.get(segment)) {
                assert(adj.fromOrientation() == -1);
                assert(adj.fromSegment() == segment);
                assert(startEdges.containsKey(adj.toSegment()));
                List<BgAdjacency> remoteAdjList;
                switch (adj.toOrientation()) {
                    case 1:
                        remoteAdjList = endEdges.get(adj.toSegment());
                        break;
                    case -1:
                        remoteAdjList = startEdges.get(adj.toSegment());
                        break;
                    default:
                        throw new IllegalArgumentException("Invalid orientation");
                }
                assert (Iterables.any(remoteAdjList, remoteAdj -> remoteAdj.edge() == adj.edge()));
            }
        }
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
