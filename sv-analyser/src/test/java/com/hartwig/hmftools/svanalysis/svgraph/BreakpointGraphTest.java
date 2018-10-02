package com.hartwig.hmftools.svanalysis.svgraph;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

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
                .build();
    }

    public EnrichedStructuralVariant breakpoint(String chr1, long pos1, int orientation1, String chr2, long pos2, int orientation2, double ploidy) {
        return ImmutableEnrichedStructuralVariant.builder()
                .id(String.format("{}:{} -> {}:{} CN:{}", chr1, pos1, chr2, pos2, ploidy))
                .type(StructuralVariantType.BND)
                .start(ImmutableEnrichedStructuralVariantLeg.builder()
                    .chromosome(chr1)
                    .position(pos1)
                    .orientation((byte)orientation1)
                    .adjustedCopyNumber(ploidy)
                    .homology("")
                    .build())
                .end(ImmutableEnrichedStructuralVariantLeg.builder()
                        .chromosome(chr2)
                        .position(pos2)
                        .orientation((byte)orientation2)
                        .adjustedCopyNumber(ploidy)
                        .homology("")
                        .build())
                .insertSequence("")
                .ploidy(1d)
                .build();
    }

    @Test
    public void BreakpointGraph_should_link_svs_at_cn_bounds() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 2),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(
                breakpoint("chr1", 10, 1, "chr1", 21, -1, 1));
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
        List<PurpleCopyNumber> cns = ImmutableList.of(
            cn("chr1", 1, 10, 2),
            cn("chr1", 11, 20, 1),
            cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(
                breakpoint("chr1", 10, 1, "chr1", 21, -1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        StructuralVariant sv = bg.simplifySimpleDeletion();
        assertEquals(svs.get(0), sv);
        assertEquals(0, bg.getAllStructuralVariants().size());
    }

    @Test
    public void simplifySimpleDeletion_should_not_remove_if_not_simple_event() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 2),
                cn("chr1", 11, 20, 1),
                cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(
                breakpoint("chr1", 10, 1, "chr1", 21, -1, 1),
                breakpoint("chr1", 11, -1, "chr1", 1, -1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        StructuralVariant sv = bg.simplifySimpleDeletion();
        assertNull(sv);
    }

    @Test
    public void simplifySimpleDuplication_should_remove_dup() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 2),
                cn("chr1", 11, 20, 3),
                cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(
                breakpoint("chr1", 11, -1, "chr1", 20, 1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        StructuralVariant sv = bg.simplifySimpleDuplications();
        assertEquals(svs.get(0), sv);
        assertEquals(0, bg.getAllStructuralVariants().size());
    }

    @Test
    public void simplifySimpleDuplication_should_not_remove_if_not_simple_event() {
        List<PurpleCopyNumber> cns = ImmutableList.of(
                cn("chr1", 1, 10, 2),
                cn("chr1", 11, 20, 3),
                cn("chr1", 21, 30, 2));
        List<EnrichedStructuralVariant> svs = ImmutableList.of(
                breakpoint("chr1", 11, -1, "chr1", 20, 1, 1),
                breakpoint("chr1", 11, -1, "chr1", 21, -1, 1));
        BreakpointGraph bg = new BreakpointGraph(cns, svs);
        StructuralVariant sv = bg.simplifySimpleDuplications();
        assertNull(sv);
    }
}