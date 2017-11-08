package com.hartwig.hmftools.common.purple.segment;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.BND;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.DEL;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.INS;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.junit.Before;
import org.junit.Test;

@Deprecated
public class PurpleSegmentFactoryOldTest {

    private GenomeRegion region1;
    private GenomeRegion region2;
    private GenomeRegion region3;
    private List<GenomeRegion> regions;

    @Before
    public void setup() {
        region1 = createRegion(1, 10000);
        region2 = createRegion(10001, 20000);
        region3 = createRegion(20001, 30000);
        regions = Lists.newArrayList(region1, region2, region3);
    }

    @Test
    public void testTooCloseToStartAndEnd() {
        final StructuralVariantPosition variant1 = createVariant(1000, StructuralVariantType.BND);
        final StructuralVariantPosition variant2 = createVariant(9000, StructuralVariantType.DEL);
        final List<StructuralVariantPosition> positions = Lists.newArrayList(variant1, variant2);
        final List<PurpleSegment> segments = PurpleSegmentFactoryOld.createSegmentsInner(regions, positions);
        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), region1, BND);
        assertPurpleSegment(segments.get(1), region2, DEL);
        assertPurpleSegment(segments.get(2), region3, NONE);
    }

    @Test
    public void testTooCloseToStart() {
        final StructuralVariantPosition variant1 = createVariant(1000, StructuralVariantType.BND);
        final StructuralVariantPosition variant2 = createVariant(7000, StructuralVariantType.INS);
        final List<StructuralVariantPosition> positions = Lists.newArrayList(variant1, variant2);
        final List<PurpleSegment> segments = PurpleSegmentFactoryOld.createSegmentsInner(regions, positions);
        assertEquals(4, segments.size());
        assertPurpleSegment(segments.get(0), 1, 6999, true, BND);
        assertPurpleSegment(segments.get(1), 7000, 10000, false, INS);
        assertPurpleSegment(segments.get(2), region2, NONE);
        assertPurpleSegment(segments.get(3), region3, NONE);
    }

    @Test
    public void testTooCloseToEnd() {
        final StructuralVariantPosition variant1 = createVariant(4000, StructuralVariantType.BND);
        final StructuralVariantPosition variant2 = createVariant(9000, StructuralVariantType.BND);
        final List<StructuralVariantPosition> positions = Lists.newArrayList(variant1, variant2);
        final List<PurpleSegment> segments = PurpleSegmentFactoryOld.createSegmentsInner(regions, positions);
        assertEquals(4, segments.size());
        assertPurpleSegment(segments.get(0), 1, 3999, true, NONE);
        assertPurpleSegment(segments.get(1), 4000, 10000, false, BND);
        assertPurpleSegment(segments.get(2), region2, BND);
        assertPurpleSegment(segments.get(3), region3, NONE);
    }

    @Test
    public void testVariantInRegion() {
        final StructuralVariantPosition variant1 = createVariant(4000, StructuralVariantType.BND);
        final StructuralVariantPosition variant2 = createVariant(7000, StructuralVariantType.BND);
        final List<StructuralVariantPosition> positions = Lists.newArrayList(variant1, variant2);
        final List<PurpleSegment> segments = PurpleSegmentFactoryOld.createSegmentsInner(regions, positions);
        assertEquals(5, segments.size());
        assertPurpleSegment(segments.get(0), 1, 3999, true, NONE);
        assertPurpleSegment(segments.get(1), 4000, 6999, false, BND);
        assertPurpleSegment(segments.get(2), 7000, 10000, false, BND);
        assertPurpleSegment(segments.get(3), region2, NONE);
        assertPurpleSegment(segments.get(4), region3, NONE);
    }

    @Test
    public void testVariantSpansRegions() {
        final StructuralVariantPosition variant1 = createVariant(3000, StructuralVariantType.BND);
        final StructuralVariantPosition variant2 = createVariant(25000, StructuralVariantType.BND);
        final List<StructuralVariantPosition> positions = Lists.newArrayList(variant1, variant2);
        final List<PurpleSegment> segments = PurpleSegmentFactoryOld.createSegmentsInner(regions, positions);
        assertEquals(5, segments.size());
        assertPurpleSegment(segments.get(0), 1, 2999, true, NONE);
        assertPurpleSegment(segments.get(1), 3000, 10000, false, BND);
        assertPurpleSegment(segments.get(2), region2, NONE);
        assertPurpleSegment(segments.get(3), 20001, 24999, true, NONE);
        assertPurpleSegment(segments.get(4), 25000, 30000, false, BND);
    }

    private void assertPurpleSegment(final PurpleSegment victim, GenomeRegion region, SegmentSupport variantSupport) {
        assertEquals(region.chromosome(), victim.chromosome());
        assertEquals(region.start(), victim.start());
        assertEquals(region.end(), victim.end());
        assertEquals(true, victim.ratioSupport());
        assertEquals(variantSupport, victim.support());
    }

    private void assertPurpleSegment(final PurpleSegment victim, long start, long end, boolean ratioSupport,
            SegmentSupport variantSupport) {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(ratioSupport, victim.ratioSupport());
        assertEquals(variantSupport, victim.support());
    }

    private static StructuralVariantPosition createVariant(long position, StructuralVariantType type) {
        return ImmutableStructuralVariantPosition.builder().chromosome("1").position(position).id("ID").orientation((byte) 1).type(type).build();
    }

    private static GenomeRegion createRegion(long start, long end) {
        return GenomeRegionFactory.create("1", start, end);
    }

}
