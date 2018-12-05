package com.hartwig.hmftools.svgraphs.simplification;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.svgraphs.BgAdjacency;
import com.hartwig.hmftools.svgraphs.BgSegment;
import com.hartwig.hmftools.svgraphs.BreakpointGraph;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

public class SimpleEventSimplifier extends Simplifier {
    private static final Logger LOGGER = LogManager.getLogger(SimpleEventSimplifier.class);
    public static final int SIMPLE_INSERTION_MAX_INSERTION_SITE_INSDEL_SIZE = 35;
    public static final int SIMPLE_INVERSION_MAX_BREAKEND_DISTANCE = 35;
    public static final int SIMPLE_INSERTION_MAX_INSERTION_LENGTH = 20000;

    public SimpleEventSimplifier(BreakpointGraph graph, SimplificationStrategy strategy) {
        super(graph, strategy);
    }

    public List<Simplification> findPotentialSimplifications() {
        List<Simplification> list = super.findPotentialSimplifications();
        list.addAll(findAssemblyLinkageSimplifications());
        list.addAll(findSimpleInversionSimplifications(SIMPLE_INVERSION_MAX_BREAKEND_DISTANCE));
        list.addAll(findInsertionSimplifications());
        list.addAll(findSimpleDuplicationSimplifications());
        list.addAll(findSimpleDeletionSimplifications());
        list.addAll(findCopyNumberLinkages());
        return list;
    }

    public void simplify(Simplification simplification) {
        assert (isValid(simplification));
        assert (strategy.shouldSimplify(simplification));
        LOGGER.debug("Simplifying {} {}", simplification.type(), simplification);
        switch (simplification.type()) {
            case SimpleDuplication:
            case SimpleIndel:
                assert (simplification.variants().size() == 1);
                Pair<BgAdjacency, BgAdjacency> adj = graph.getAdjacencies(simplification.variants().get(0));
                BgAdjacency left = adj.getLeft();
                BgAdjacency right = adj.getRight();
                List<BgSegment> segments = graph.getSegmentsBetween(left, right);
                assert (segments.size() == 1);
                graph.removeEdge(left, right);
                for (BgSegment segment : segments) {
                    // removing a deletion adds to the CN
                    // removing a duplication reduces the CN
                    segment.adjustCopyNumber(simplification.copyNumber() * (simplification.type() == SimplificationType.SimpleIndel ? 1 : -1));
                }
                return;
            case TranslocationInsertion:
                simplifyTranslocationInsertion(simplification);
                return;
            case Chain:
                simplifyChain(simplification);
                return;
            case SimpleInversion:
                simplifyInversion(simplification);
                return;
            default:
                super.simplify(simplification);
        }
    }

    private void simplifyInversion(Simplification simplification) {
        BgAdjacency left1 = simplification.adjacencies().get(0);
        BgAdjacency left2 = simplification.adjacencies().get(1);
        BgAdjacency right1 = simplification.adjacencies().get(2);
        BgAdjacency right2 = simplification.adjacencies().get(3);
        assert(!left1.isReference());
        assert(left1.fromOrientation() == left1.toOrientation());
        assert(left2.fromOrientation() == left2.toOrientation());
        assert(left1.fromOrientation() != left2.fromOrientation());
        assert(left1.edge().sv() != left2.edge().sv());
        // Adjust CN of duplicated/lost sections
        // left-most SV has a + orientation -> boundary sequence lost
        // left-most SV has a - orientation -> boundary sequence duplicated
        List<BgSegment> leftSegments = graph.getSegmentsBetween(left1, left2);
        for (BgSegment seg : leftSegments) {
            seg.adjustCopyNumber(simplification.copyNumber() * (left1.fromOrientation() == 1 ? 1 : -1));
        }
        List<BgSegment> rightSegments = graph.getSegmentsBetween(right1, right2);
        for (BgSegment seg : rightSegments) {
            seg.adjustCopyNumber(simplification.copyNumber() * (right1.fromOrientation() == 1 ? 1 : -1));
        }
        graph.removeEdge(left1);
        graph.removeEdge(left2);
    }

    private void simplifyChain(Simplification simplification) {
        assert (simplification.adjacencies().size() == 2);
        graph.mergeAdjacentSvs(simplification.copyNumber(), simplification.adjacencies().get(0), simplification.adjacencies().get(1));
    }

    private void simplifyTranslocationInsertion(Simplification simplification) {
        assert (simplification.variants().size() == 2 || simplification.variants().size() == 3);
        assert (simplification.variants().stream().distinct().count() == simplification.variants().size());
        BgAdjacency left = simplification.adjacencies().get(0);
        BgAdjacency right = simplification.adjacencies().get(1);
        List<BgSegment> segments = graph.getSegmentsBetween(left, right);
        for (BgSegment segment : segments) {
            // -ve on the left side indicates that the sequence is duplicated
            segment.adjustCopyNumber(simplification.copyNumber() * (left.fromOrientation() == -1 ? -1 : 1));
        }
        BgAdjacency remoteA = graph.getPartner(left);
        BgAdjacency remoteB = graph.getPartner(right);
        if (remoteA.fromSegment() == graph.getUnplacedSegment()) {
            remoteA = simplification.adjacencies().size() == 3 ? simplification.adjacencies().get(2) : null;
        } else if (remoteB.fromSegment() == graph.getUnplacedSegment()) {
            remoteB = simplification.adjacencies().size() == 3 && remoteA != simplification.adjacencies().get(2) ? simplification.adjacencies().get(2) : null;
        }
        if (remoteA != null && remoteB != null) {
            List<BgSegment> insertedSegment = graph.getSegmentsBetween(remoteA, remoteB);
            for (BgSegment segment : insertedSegment) {
                // -ve on the left side indicates that the sequence is duplicated
                segment.adjustCopyNumber(-simplification.copyNumber());
            }
        } else {
            // TODO: what do we do if we only have 1 side of the insertion?
        }
        for (BgAdjacency adj : simplification.adjacencies()) {
            graph.removeEdge(adj);
        }
    }

    public List<Simplification> findSimpleDuplicationSimplifications() {
        List<Simplification> result = new ArrayList<>();
        for (BgSegment segment : graph.getAllSegments()) {
            List<BgAdjacency> prev = graph.getOutgoing(segment, -1);
            List<BgAdjacency> next = graph.getOutgoing(segment, 1);
            // simple duplication check
            if (graph.svCount(prev) == 1 && graph.svCount(next) == 1) {
                BgAdjacency prevSv = graph.firstSv(prev);
                BgAdjacency nextSv = graph.firstSv(next);
                if (prevSv.edge() == nextSv.edge()) {
                    BgSegment prevSegment = graph.prevReferenceSegment(segment);
                    BgSegment nextSegment = graph.nextReferenceSegment(segment);
                    if (graph.svCount(graph.getOutgoing(prevSegment, 1)) == 0 && graph.svCount(graph.getOutgoing(nextSegment, -1)) == 0) {
                        EnrichedStructuralVariant sv = prevSv.edge().sv();
                        double copyNumber = segment.copyNumber() - (prevSegment.copyNumber() + nextSegment.copyNumber()) / 2;
                        Simplification duplication = ImmutableSimplificationImpl.builder()
                                .type(SimplificationType.SimpleDuplication)
                                .copyNumber(copyNumber)
                                .addVariants(sv)
                                .consistency(graph.getConsistencySet(copyNumber, nextSv))
                                .addAdjacencies(prevSv, nextSv)
                                .build();
                        if (strategy.shouldSimplify(duplication)) {
                            result.add(duplication);
                        }
                    }
                }
            }
        }
        return result;
    }

    public List<Simplification> findSimpleDeletionSimplifications() {
        List<Simplification> result = new ArrayList<>();
        for (BgSegment segment : graph.getAllSegments()) {
            BgSegment prevSegment = graph.prevReferenceSegment(segment);
            BgSegment nextSegment = graph.nextReferenceSegment(segment);
            if (prevSegment != null && nextSegment != null && graph.svCount(graph.getOutgoing(segment, -1)) == 0
                    && graph.svCount(graph.getOutgoing(segment, 1)) == 0 && graph.svCount(graph.getOutgoing(prevSegment, 1)) == 1
                    && graph.svCount(graph.getOutgoing(nextSegment, -1)) == 1) {
                BgAdjacency leftSv = graph.firstSv(graph.getOutgoing(prevSegment, 1));
                BgAdjacency rightSv = graph.firstSv(graph.getOutgoing(nextSegment, -1));
                if (leftSv.edge() == rightSv.edge()) {
                    // we have a deletion
                    double copyNumber = (prevSegment.copyNumber() + nextSegment.copyNumber()) / 2 - segment.copyNumber();
                    Simplification deletion = ImmutableSimplificationImpl.builder()
                            .type(SimplificationType.SimpleIndel)
                            .copyNumber(copyNumber)
                            .addVariants(leftSv.edge().sv())
                            .consistency(graph.getConsistencySet(copyNumber, leftSv))
                            .addAdjacencies(leftSv, rightSv)
                            .build();
                    if (strategy.shouldSimplify(deletion)) {
                        result.add(deletion);
                    }
                }
            }
        }
        return result;
    }

    public List<Simplification> findAssemblyLinkageSimplifications() {
        return findPreLinkedVariants("asm");
    }

    public List<Simplification> findPreLinkedVariants(String linkedByPrefix) {
        return findAdjacentBreakendLinks(graph.getLinkedBreakendPairs(linkedByPrefix));
    }

    public List<Simplification> findAdjacentBreakendLinks(Collection<Pair<BgAdjacency, BgAdjacency>> links) {
        List<Simplification> result = new ArrayList<>();
        for (Pair<BgAdjacency, BgAdjacency> link : links) {
            BgAdjacency left = link.getLeft();
            BgAdjacency right = link.getRight();
            double copyNumber = (left.edge().ploidy() + right.edge().ploidy()) / 2;
            Simplification chain = ImmutableSimplificationImpl.builder()
                    .type(SimplificationType.Chain)
                    .copyNumber(copyNumber)
                    .addVariants(left.edge().sv(), right.edge().sv())
                    .consistency(graph.getConsistencySet(copyNumber, left, right))
                    .addAdjacencies(left, right)
                    .build();
            if (strategy.shouldSimplify(chain)) {
                result.add(chain);
            }
        }
        return result;
    }

    public List<Simplification> findCopyNumberLinkages() {
        List<Simplification> result = new ArrayList<>();
        for (BgSegment segment : graph.getAllSegments()) {
            if (segment != graph.getUnplacedSegment() && graph.svCount(graph.getOutgoing(segment, -1)) == 1) {
                BgSegment prev = graph.prevReferenceSegment(segment);
                if (prev != null && graph.svCount(graph.getOutgoing(prev, 1)) == 0) {
                    BgAdjacency adj = graph.firstSv(graph.getOutgoing(segment, -1));
                    List<BgAdjacency> foldBackDoublingPartners = graph.nextFoldbackDoublingCandidates(adj);
                    if (foldBackDoublingPartners.size() > 0) {
                        // TODO: process fold-back inversions
                        continue;
                        // need to check that the copy number between the two foldback breakend is consistent with this event
                        // being directly linked to the foldback inversion
                        //if (strategy.couldBeFoldBackLinked(adj.edge().sv(), getSegmentsBetween(foldback, getPartner(foldback)))) {
                        //}
                    }
                    List<BgAdjacency> partners = graph.nextBreakpointCandidates(adj, null);
                    if (partners.size() == 1 && partners.get(0) != null) {
                        if (strategy.couldBeFoldBackLinked(adj.edge().sv(), partners.get(0).edge().sv())) {
                            // Don't simplify across fold-back linkage as we halve/double our copy number
                            // when we traverse back across a fold-back inversion. If we simplify the event
                            // then we get our copyNumber incorrect
                            LOGGER.debug("Not simplifying to potential fold-back inversion {} to {}",
                                    adj.edge().sv(),
                                    partners.get(0).edge().sv());
                        } else {
                            BgAdjacency partner = partners.get(0);
                            if (partner.edge().sv() != adj.edge().sv()) {
                                double copyNumber = (adj.edge().ploidy() + partner.edge().ploidy()) / 2;
                                Simplification chain = ImmutableSimplificationImpl.builder()
                                        .type(SimplificationType.Chain)
                                        .copyNumber(copyNumber)
                                        .addVariants(adj.edge().sv(), partner.edge().sv())
                                        .consistency(graph.getConsistencySet(copyNumber, graph.getPartner(adj), partner))
                                        .addAdjacencies(partner.fromOrientation() == -1 ? partner : graph.getPartner(adj),
                                                partner.fromOrientation() == -1 ? graph.getPartner(adj) : partner)
                                        .build();
                                if (strategy.shouldSimplify(chain)) {
                                    result.add(chain);
                                }
                            }
                        }
                    }
                }
            }
        }
        return result;
    }
    public List<Simplification> findInsertionSimplifications() {
        List<Simplification> result = new ArrayList<>();
        for (Pair<BgAdjacency, BgAdjacency> pair : graph.findPotentialSimpleInsertionSites(SIMPLE_INSERTION_MAX_INSERTION_SITE_INSDEL_SIZE, SIMPLE_INSERTION_MAX_INSERTION_LENGTH)) {
            BgAdjacency left = pair.getLeft();
            BgAdjacency right = pair.getRight();
            Simplification simplification = createInsertionSimplification(left, right);
            if (strategy.shouldSimplify(simplification)) {
                result.add(simplification);
            }
        }
        return result;
    }
    public List<Simplification> findSimpleInversionSimplifications(long simpleInversionMaxBreakendDistance) {
        List<Simplification> result = new ArrayList<>();
        List<BgAdjacency> adjList = graph.getAllAdjacencies()
                .filter(adj -> !adj.isReference())
                .sorted(graph.ByFromGenomicPosition)
                .collect(Collectors.toList());
        for (int i = 0; i < adjList.size() - 3; i++) {
            BgAdjacency a = adjList.get(i);
            BgAdjacency b = adjList.get(i + 1);
            BgAdjacency c = adjList.get(i + 2);
            BgAdjacency d = adjList.get(i + 3);
            if (a.fromSegment().chromosome().equals(d.fromSegment().chromosome()) &&
                    a.fromOrientation() == a.toOrientation() &&
                    b.fromOrientation() == b.toOrientation() &&
                    a.fromOrientation() != b.toOrientation() &&
                    ((a.edge() == c.edge() && b.edge() == d.edge()) || a.edge() == d.edge() && b.edge() == c.edge())) {
                long leftImperfectionLength = graph.getSegmentsBetween(a, b).stream().mapToLong(s -> s.length()).sum();
                long rightImperfectionLength  = graph.getSegmentsBetween(a, b).stream().mapToLong(s -> s.length()).sum();
                if (leftImperfectionLength <= simpleInversionMaxBreakendDistance &&  rightImperfectionLength <= simpleInversionMaxBreakendDistance) {
                    double copyNumber = (a.edge().ploidy() + b.edge().ploidy()) / 2;
                    ImmutableSimplificationImpl.Builder builder = ImmutableSimplificationImpl.builder()
                            .type(SimplificationType.SimpleInversion)
                            .copyNumber(copyNumber)
                            .addVariants(a.edge().sv(), b.edge().sv())
                            .consistency(graph.getConsistencySet(copyNumber, a, b))
                            .addAdjacencies(a, b, c, d);
                    Simplification simplification = builder.build();
                    if (strategy.shouldSimplify(simplification)) {
                        result.add(simplification);
                    }
                }
            }
        }
        return result;
    }

    private Simplification createInsertionSimplification(BgAdjacency left, BgAdjacency right) {
        double copyNumber = (left.edge().ploidy() + right.edge().ploidy()) / 2;
        ImmutableSimplificationImpl.Builder builder = ImmutableSimplificationImpl.builder()
                .type(SimplificationType.TranslocationInsertion)
                .copyNumber(copyNumber)
                .addVariants(left.edge().sv(), right.edge().sv())
                // do we care about the consistency of the flanking region?
                .consistency(graph.getConsistencySet(copyNumber, left, right))
                .addAdjacencies(left, right);
        if ((left.toSegment() == graph.getUnplacedSegment() && right.toSegment() != graph.getUnplacedSegment()) ||
                (left.toSegment() != graph.getUnplacedSegment() && right.toSegment() == graph.getUnplacedSegment())) {
            // attempt to associate any dangling single breakend at the site of the insertion
            // to complete our templated insertion
            BgAdjacency bp = left.toSegment() == graph.getUnplacedSegment() ? right : left;
            BgAdjacency bpInsFlank = graph.getPartner(left.toSegment() == graph.getUnplacedSegment() ? right : left);
            List<BgAdjacency> potentialEnds = graph.nextBreakpointCandidates(bp, null).stream()
                    .filter(x -> x != null
                            && x.toSegment() == graph.getUnplacedSegment()
                            && graph.segmentLength(bpInsFlank, x) != null
                            && graph.segmentLength(bpInsFlank, x) <= SIMPLE_INSERTION_MAX_INSERTION_LENGTH)
                    .collect(Collectors.toList());
            if (potentialEnds.size() == 1) {
                builder = builder.addVariants(potentialEnds.get(0).edge().sv())
                        .addAdjacencies(potentialEnds.get(0));
            }
        }
        Simplification simplification = builder.build();
        return simplification;
    }
}
