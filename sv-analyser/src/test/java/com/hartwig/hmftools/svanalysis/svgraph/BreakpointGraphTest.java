package com.hartwig.hmftools.svanalysis.svgraph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Stream;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

public class BreakpointGraphTest {
    public PurpleCopyNumber cn(String chr, long start, long end, double ploidy) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chr)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(ploidy)
                .method(CopyNumberMethod.BAF_WEIGHTED)
                .bafCount(0)
                .averageActualBAF(0)
                .averageObservedBAF(0)
                .depthWindowCount(1)
                .segmentStartSupport(SegmentSupport.UNKNOWN)
                .segmentEndSupport(SegmentSupport.UNKNOWN)
                .gcContent(0)
                .minStart(start)
                .maxStart(start)
                .build();
    }

    public EnrichedStructuralVariant breakpoint(String chr1, long pos1, int orientation1, String chr2, long pos2, int orientation2,
                                                double ploidy) {
        return breakpoint(String.format("{}:{} -> {}:{} CN:{}", chr1, pos1, chr2, pos2, ploidy),
                chr1,
                pos1,
                orientation1,
                chr2,
                pos2,
                orientation2,
                ploidy);
    }

    public EnrichedStructuralVariant breakpoint(String id, String chr1, long pos1, int orientation1, String chr2, long pos2,
                                                int orientation2, double ploidy) {
        return ImmutableEnrichedStructuralVariant.builder()
                .id(id)
                .type(StructuralVariantType.BND)
                .start(ImmutableEnrichedStructuralVariantLeg.builder()
                        .chromosome(chr1)
                        .position(pos1)
                        .orientation((byte) orientation1)
                        .adjustedCopyNumber(ploidy)
                        .homology("")
                        .build())
                .end(ImmutableEnrichedStructuralVariantLeg.builder()
                        .chromosome(chr2)
                        .position(pos2)
                        .orientation((byte) orientation2)
                        .adjustedCopyNumber(ploidy)
                        .homology("")
                        .build())
                .insertSequence("")
                .ploidy(ploidy)
                .qualityScore(1000.0)
                .recovered(false)
                .startLinkedBy("")
                .endLinkedBy("")
                .filter("")
                .imprecise(false)
                .insertSequence("")
                .build();
    }
    public EnrichedStructuralVariant breakend(String id, String chr1, long pos1, int orientation1, double ploidy) {
        return ImmutableEnrichedStructuralVariant.builder()
                .id(id)
                .type(StructuralVariantType.BND)
                .start(ImmutableEnrichedStructuralVariantLeg.builder()
                        .chromosome(chr1)
                        .position(pos1)
                        .orientation((byte) orientation1)
                        .adjustedCopyNumber(ploidy)
                        .homology("")
                        .build())
                .insertSequence("")
                .ploidy(ploidy)
                .qualityScore(1000.0)
                .recovered(false)
                .startLinkedBy("")
                .endLinkedBy("")
                .filter("")
                .imprecise(false)
                .insertSequence("GATTACA")
                .build();
    }

    @Test
    public void BreakpointGraph_should_link_svs_at_cn_bounds() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 2), cn("chr1", 11, 20, 1), cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(breakpoint("chr1", 10, 1, "chr1", 21, -1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        assertEquals(4, bg.getAllSegments().size());
        assertEquals(3, bg.getAllEdges().size());
        assertEquals(1, bg.getAllStructuralVariants().size());
        assertEquals(0, bg.getOutgoing(bg.getUnplacedSegment(), 1).size());
        assertEquals(0, bg.getOutgoing(bg.getUnplacedSegment(), -1).size());
        assertEquals(2, bg.getOutgoing(bg.getSegment(cns.get(0)), 1).size());
        assertEquals(2, bg.getOutgoing(bg.getSegment(cns.get(2)), -1).size());
    }

    @Test
    public void simplifySimpleDeletion_should_remove_del() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 2), cn("chr1", 11, 20, 1), cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(breakpoint("chr1", 10, 1, "chr1", 21, -1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        Simplification sv = bg.simplify().get(0);
        assertEquals(svs.get(0), sv.variants().get(0));
        assertEquals(0, bg.getAllStructuralVariants().size());
    }

    @Test
    public void simplifySimpleDeletion_should_not_remove_if_not_simple_event() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 2), cn("chr1", 11, 20, 1), cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs =
                ImmutableList.of(breakpoint("chr1", 10, 1, "chr1", 21, -1, 1), breakpoint("chr1", 11, -1, "chr1", 1, -1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        List<Simplification> results = bg.simplify();
        assertEquals(0, results.size());
    }

    @Test
    public void simplifySimpleDuplication_should_remove_dup() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 2), cn("chr1", 11, 20, 3), cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(breakpoint("chr1", 11, -1, "chr1", 20, 1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        Simplification sv = bg.simplify().get(0);
        assertEquals(svs.get(0), sv.variants().get(0));
        assertEquals(0, bg.getAllStructuralVariants().size());
    }

    @Test
    public void simplifySimpleDuplication_should_not_remove_if_not_simple_event() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 2), cn("chr1", 11, 20, 3), cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs =
                ImmutableList.of(breakpoint("chr1", 11, -1, "chr1", 20, 1, 1), breakpoint("chr1", 11, -1, "chr1", 21, -1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        List<Simplification> results = bg.simplify();
        assertEquals(0, results.size());
    }

    @Test
    public void simplifyAssemblyLinks() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 2),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 2),
                cn("chr1", 31, 40, 1),
                cn("chr1", 41, 50, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(ImmutableEnrichedStructuralVariant.builder()
                        .from(breakpoint("chr1", 10, 1, "chr1", 21, -1, 1))
                        .startLinkedBy("asm_left")
                        .endLinkedBy("asm001")
                        .build(),
                ImmutableEnrichedStructuralVariant.builder()
                        .from(breakpoint("chr1", 30, 1, "chr1", 41, -1, 1))
                        .startLinkedBy("asm001")
                        .endLinkedBy("asm_right")
                        .build());
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        List<Simplification> sv = bg.simplify();
        assertEquals(1, sv.size());
        assertEquals(2, sv.get(0).variants().size());
        assertEquals(1, bg.getAllStructuralVariants().size());
        assertEquals(10, bg.getAllStructuralVariants().get(0).start().position());
        assertEquals(1, bg.getAllStructuralVariants().get(0).start().orientation());
        assertEquals(41, bg.getAllStructuralVariants().get(0).end().position());
        assertEquals(-1, bg.getAllStructuralVariants().get(0).end().orientation());
        assertEquals("asm_left", bg.getAllStructuralVariants().get(0).startLinkedBy());
        assertEquals("asm_right", bg.getAllStructuralVariants().get(0).endLinkedBy());
    }

    @Test
    public void candidatesNextBreakpoint_should_find_potential_partners() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 2),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 2),
                cn("chr1", 31, 40, 1),
                cn("chr1", 41, 50, 2),
                cn("other", 1, 1, 1));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(breakpoint("sv1", "other", 1, -1, "chr1", 1, -1, 1),
                breakpoint("sv2", "other", 1, 1, "chr1", 10, 1, 1),
                breakpoint("sv3", "other", 1, 1, "chr1", 20, 1, 1),
                breakpoint("sv4", "other", 1, 1, "chr1", 30, 1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        BgAdjacency ba = bg.getAllAdjacencies()
                .filter(adj -> adj.fromSegment().chromosome().equals("other") && adj.edge().sv() != null && adj.edge()
                        .sv()
                        .id()
                        .equals("sv1"))
                .findFirst()
                .orElse(null);
        assertNotNull(ba);
        List<BgAdjacency> result = bg.nextBreakpointCandidates(ba);
        assertEquals(4, result.size());
    }

    @Test
    public void candidatesNextBreakpoint_should_not_match_if_CN_too_different() {
        List<PurpleCopyNumber> cns = ImmutableList.of(cn("chr1", 1, 10, 10),
                cn("chr1", 11, 20, 10),
                cn("chr1", 21, 30, 1),
                cn("chr1", 31, 40, 10),
                cn("chr1", 41, 50, 1),
                cn("other", 1, 1, 1));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(breakpoint("sv1", "other", 1, -1, "chr1", 1, -1, 5),
                breakpoint("sv2", "other", 1, 1, "chr1", 10, 1, 1),
                breakpoint("sv3", "other", 1, 1, "chr1", 20, 1, 5),
                breakpoint("sv4", "other", 1, 1, "chr1", 30, 1, 5));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        BgAdjacency ba = bg.getAllAdjacencies()
                .filter(adj -> adj.fromSegment().chromosome().equals("other") && adj.edge().sv() != null && adj.edge()
                        .sv()
                        .id()
                        .equals("sv1"))
                .findFirst()
                .orElse(null);
        assertNotNull(ba);
        List<BgAdjacency> result = bg.nextBreakpointCandidates(ba);
        // sv3 is the only candidate left
        // sv2 has too low copyNumber
        // can't traverse to either sv4 or end of chromosome due to low CN
        assertEquals(1, result.size());
        assertEquals("sv3", result.get(0).edge().sv().id());
    }

    @Test
    public void nextFoldbackDoublingCandidates() {
        List<PurpleCopyNumber> cns =
                ImmutableList.of(cn("chr1", 1, 10, 11), cn("chr1", 11, 20, 11), cn("chr1", 21, 30, 6), cn("other", 1, 1, 1));
        List<EnrichedStructuralVariant> svs =
                ImmutableList.of(breakpoint("sv1", "other", 1, -1, "chr1", 1, -1, 10), breakpoint("sv2", "chr1", 30, 1, "chr1", 20, 1, 5));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        BgAdjacency ba = bg.getAllAdjacencies()
                .filter(adj -> adj.fromSegment().chromosome().equals("other") && adj.edge().sv() != null && adj.edge()
                        .sv()
                        .id()
                        .equals("sv1"))
                .findFirst()
                .orElse(null);
        assertNotNull(ba);
        List<BgAdjacency> result = bg.nextFoldbackDoublingCandidates(ba);
        // can follow a fold-back inversion from either direction
        assertEquals(2, result.size());
    }

    @Test
    public void findPotentialSimpleInsertionSites_breakend() {
        //
        //  chr1   AAA BBB CCC
        //  chr2   XXX YYY ZZZ
        //
        // Insertion of BBB at XY/YZ
        //
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 1),
                cn("chr2", 1, 100, 1),
                cn("chr2", 101, 200, 2),
                cn("chr2", 201, 300, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 11, -1, "chr2", 101, -1, 1),
                breakpoint("right", "chr1", 20, 1, "chr2", 200, 1, 1)));
        List<Pair<BgAdjacency, BgAdjacency>> sites = bg.findPotentialSimpleInsertionSites(20, 1000);
        assertEquals(1, sites.size());
        assertEquals(-1, sites.get(0).getLeft().fromOrientation());
        assertEquals("chr1", sites.get(0).getLeft().fromSegment().chromosome());
        assertEquals(1, sites.get(0).getRight().fromOrientation());
        assertEquals("chr1", sites.get(0).getRight().fromSegment().chromosome());

        // inverted insertion should also be fine
        bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 20, 1, "chr2", 101, -1, 1),
                breakpoint("right", "chr1", 11, -1, "chr2", 200, 1, 1)));
        sites = bg.findPotentialSimpleInsertionSites(20, 1000);
        assertEquals(1, sites.size());
        assertEquals(-1, sites.get(0).getLeft().fromOrientation());
        assertEquals("chr1", sites.get(0).getLeft().fromSegment().chromosome());
        assertEquals(1, sites.get(0).getRight().fromOrientation());
        assertEquals("chr1", sites.get(0).getRight().fromSegment().chromosome());
    }

    @Test
    public void findPotentialSimpleInsertionSites_should_respect_max_size_bounds() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 1),
                cn("chr2", 1, 100, 1),
                cn("chr2", 101, 200, 2),
                cn("chr2", 201, 300, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 11, -1, "chr2", 101, -1, 1),
                breakpoint("right", "chr1", 20, 1, "chr2", 200, 1, 1)));
        assertEquals(1, bg.findPotentialSimpleInsertionSites(11, 101).size());
        assertEquals(1, bg.findPotentialSimpleInsertionSites(11, 100).size());
        assertEquals(1, bg.findPotentialSimpleInsertionSites(10, 101).size());
        assertEquals(1, bg.findPotentialSimpleInsertionSites(10, 100).size());
        assertEquals(0, bg.findPotentialSimpleInsertionSites(9, 100).size());
        assertEquals(0, bg.findPotentialSimpleInsertionSites(10, 99).size());
    }
    @Test
    public void findPotentialSimpleInsertionSites_breakpoint() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 1),
                cn("chr2", 1, 100, 1),
                cn("chr2", 101, 200, 2),
                cn("chr2", 201, 300, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 11, -1, "chr2", 101, -1, 1),
                breakend("right", "chr1", 20, 1, 1)));
        List<Pair<BgAdjacency, BgAdjacency>> sites = bg.findPotentialSimpleInsertionSites(1000, 100);
        assertEquals(1, sites.size());
        assertEquals(-1, sites.get(0).getLeft().fromOrientation());
        assertEquals("chr1", sites.get(0).getLeft().fromSegment().chromosome());
        assertEquals(1, sites.get(0).getRight().fromOrientation());
        assertEquals("chr1", sites.get(0).getRight().fromSegment().chromosome());
    }
    @Test
    public void findPotentialSimpleInsertionSites_should_not_self_match() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 2),
                cn("chr1", 21, 30, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("dup", "chr1", 11, -1, "chr1", 20, 1, 1)));
        List<Pair<BgAdjacency, BgAdjacency>> sites = bg.findPotentialSimpleInsertionSites(1000, 100);
        assertEquals(0, sites.size());
    }
    @Test
    public void simplifySimpleInsertionTranslocation_breakpoint() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 3),
                cn("chr1", 21, 30, 1),
                cn("chr2", 1, 100, 1),
                cn("chr2", 101, 200, 3),
                cn("chr2", 201, 300, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 11, -1, "chr2", 101, -1, 2),
                breakpoint("right", "chr1", 20, 1, "chr2", 200, 1, 2)));
        List<Simplification> simplifications = bg.simplify();
        assertEquals(1, simplifications.size());
        assertEquals(SimplificationType.TranslocationInsertion, simplifications.get(0).type());

        bg.mergeReferenceSegments(true);
        assertEquals(3, bg.getAllSegments().size());
    }
    @Test
    public void simplifySimpleInsertionTranslocation_should_rescue_breakend() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 2),
                cn("chr1", 21, 30, 1),
                cn("chr2", 1, 100, 1),
                cn("chr2", 101, 200, 2),
                cn("chr2", 201, 300, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 11, -1, "chr2", 101, -1, 1),
                breakend("right", "chr1", 20, 1, 1),
                breakend("remote", "chr2", 200, 1, 1)));
        List<Simplification> simplifications = bg.simplify();
        bg.mergeReferenceSegments(true);
        assertEquals(1, simplifications.size());
        assertEquals(3, simplifications.get(0).variants().size());
        assertTrue(bg.getAllSegments().stream()
                .filter(x -> x != bg.getUnplacedSegment())
                .allMatch(x -> x.copyNumber() == 1));
    }
    private BgAdjacency adj(BreakpointGraph bg, String id, String chr, Long position, Integer orientation) {
        Stream<BgAdjacency> s = bg.getAllAdjacencies();
        if (id != null) {
            s = s.filter(adj -> !adj.isReference() && adj.edge().sv().id().equals(id));
        }
        if (chr != null) {
            s = s.filter(adj -> adj.fromSegment().chromosome().equals(chr));
        }
        if (position != null) {
            s = s.filter(adj -> adj.fromSegment().startPosition() == position || adj.fromSegment().endPosition() == position);
        }
        if (orientation != null) {
            s = s.filter(adj -> adj.fromOrientation() == orientation);
        }
        BgAdjacency adj = s.findFirst().orElse(null);
        return adj;
    }
    @Test
    public void segmentLength() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 2),
                cn("chr1", 21, 30, 1),
                cn("chr2", 1, 100, 1),
                cn("chr2", 101, 200, 2),
                cn("chr2", 201, 300, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("sv", "chr1", 11, -1, "chr2", 101, -1, 1),
                breakend("be", "chr1", 20, 1, 1),
                breakpoint("bp", "chr1", 20, 1, "chr2", 101, -1, 1),
                breakend("wrongChr", "chr2", 200, 1, 1),
                breakend("wrongDir", "chr1", 21, -1, 1),
                breakend("wrongPos", "chr1", 10, 1, 1)));
        BgAdjacency a = adj(bg, "sv", "chr1", null, null);
        assertEquals((Long)10L, bg.segmentLength(a, adj(bg, "be", "chr1", null, null)));
        assertEquals((Long)10L, bg.segmentLength(a, adj(bg, "bp", "chr1", null, null)));
        assertEquals(null, bg.segmentLength(a, adj(bg, "wrongChr", null, null, null)));
        assertEquals(null, bg.segmentLength(a, adj(bg, "wrongDir", null, null, null)));
        assertEquals(null, bg.segmentLength(a, adj(bg, "wrongPos", null, null, null)));
    }
    @Test
    public void simplifyInversion_gain_loss() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 2),
                cn("chr1", 21, 30, 1),
                cn("chr1", 31, 40, 0),
                cn("chr1", 41, 50, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 11, -1, "chr1", 31, -1, 1),
                breakpoint("right", "chr1", 20, 1, "chr1", 40, 1, 1)));
        List<Simplification> simplifications = bg.findSimpleInversionSimplifications(9);
        assertEquals(0, simplifications.size());
        simplifications = bg.findSimpleInversionSimplifications(10);
        assertEquals(1, simplifications.size());
        simplifications = bg.simplify();
        bg.mergeReferenceSegments(false);
        assertEquals(1, simplifications.size());
        assertEquals(SimplificationType.SimpleInversion, simplifications.get(0).type());
        assertEquals(2, bg.getAllSegments().size());
    }
    @Test
    public void simplifyInversion_perfect_breaks() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 11, -1, "chr1", 21, -1, 1),
                breakpoint("right", "chr1", 10, 1, "chr1", 20, 1, 1)));
        List<Simplification> simplifications = bg.findSimpleInversionSimplifications(10);
        assertEquals(1, simplifications.size());
        simplifications = bg.simplify();
        bg.mergeReferenceSegments(false);
        assertEquals(1, simplifications.size());
        assertEquals(SimplificationType.SimpleInversion, simplifications.get(0).type());
        assertEquals(2, bg.getAllSegments().size());
    }
    @Test
    public void simplifyInversion_loss_loss() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 0),
                cn("chr1", 21, 30, 1),
                cn("chr1", 31, 40, 0),
                cn("chr1", 41, 50, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                breakpoint("left", "chr1", 21, -1, "chr1", 31, -1, 1),
                breakpoint("right", "chr1", 10, 1, "chr1", 40, 1, 1)));
        List<Simplification> simplifications = bg.findSimpleInversionSimplifications(10);
        assertEquals(1, simplifications.size());
        simplifications = bg.simplify();
        bg.mergeReferenceSegments(false);
        assertEquals(1, simplifications.size());
        assertEquals(SimplificationType.SimpleInversion, simplifications.get(0).type());
        assertEquals(2, bg.getAllSegments().size());
    }
    @Test
    public void mergeAdjacentSvs_should_retain_links() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 1),
                cn("chr1", 11, 20, 0),
                cn("chr1", 21, 30, 1),
                cn("chr1", 31, 40, 0),
                cn("chr1", 41, 50, 1),
                cn("chr1", 51, 60, 0),
                cn("chr1", 61, 70, 1));

        BreakpointGraph bg = new BreakpointGraph(cns, ImmutableList.of(
                ImmutableEnrichedStructuralVariant.builder().from(
                    breakpoint("merge_left", "chr1", 10, 1, "chr1", 21, -1, 1))
                    .endLinkedBy("asm1").build(),
                ImmutableEnrichedStructuralVariant.builder().from(
                        breakpoint("merge_right", "chr1", 30, 1, "chr1", 41, -1, 1))
                        .startLinkedBy("asm1").endLinkedBy("asm2").build(),
                ImmutableEnrichedStructuralVariant.builder().from(
                        breakpoint("chain", "chr1", 50, 1, "chr1", 61, -1, 1))
                        .startLinkedBy("asm2").build()));
        bg.simplify();
    }
}
