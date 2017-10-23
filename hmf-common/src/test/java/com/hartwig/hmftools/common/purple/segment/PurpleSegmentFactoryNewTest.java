package com.hartwig.hmftools.common.purple.segment;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.ImmutableChromosomeLength;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleSegmentFactoryNewTest {
    private static final ChromosomeLength CHROM = ImmutableChromosomeLength.builder().chromosome("chromosome").position(10_000_000).build();

    @Test
    public void testEmpty() {
        final List<PurpleSegment> segments = PurpleSegmentFactoryNew.create(CHROM, Collections.emptyList(), Collections.emptyList());
        assertEquals(1, segments.size());
        assertPurpleSegment(segments.get(0), 1, CHROM.position(), true, StructuralVariantSupport.NONE);
    }

    @Test
    public void testSingleSV() {
        final List<StructuralVariantCluster> clusters = Lists.newArrayList(builder(18881, 17001, 20000).build());
        final List<PurpleSegment> segments = PurpleSegmentFactoryNew.create(CHROM, clusters, Collections.emptyList());
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, StructuralVariantSupport.NONE);
        assertPurpleSegment(segments.get(1), 18881, CHROM.position(), false, StructuralVariantSupport.BND);
    }

    @Test
    public void testSingleSVWithRatioSupportAtStart() {
        final List<StructuralVariantCluster> clusters = Lists.newArrayList(builder(18881, 17001, 20000).build());
        final List<GenomePosition> ratios = Lists.newArrayList(ratio(17050));
        final List<PurpleSegment> segments = PurpleSegmentFactoryNew.create(CHROM, clusters, ratios);
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, StructuralVariantSupport.NONE);
        assertPurpleSegment(segments.get(1), 18881, CHROM.position(), true, StructuralVariantSupport.BND);
    }

    @Test
    public void testSingleSVWithRatioSupportAtEnd() {
        final List<StructuralVariantCluster> clusters = Lists.newArrayList(builder(18881, 17001, 20000).build());
        final List<GenomePosition> ratios = Lists.newArrayList(ratio(19050));
        final List<PurpleSegment> segments = PurpleSegmentFactoryNew.create(CHROM, clusters, ratios);
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, StructuralVariantSupport.NONE);
        assertPurpleSegment(segments.get(1), 18881, CHROM.position(), true, StructuralVariantSupport.BND);
    }

    private static void assertPurpleSegment(final PurpleSegment victim, long start, long end, boolean ratioSupport,
            StructuralVariantSupport variantSupport) {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(ratioSupport, victim.ratioSupport());
        assertEquals(variantSupport, victim.structuralVariantSupport());
    }

    private static ImmutableStructuralVariantCluster.Builder builder(long position, long start, long end) {
        return ImmutableStructuralVariantCluster.builder()
                .chromosome(CHROM.chromosome())
                .start(start)
                .end(end)
                .addVariants(createSVPosition(position));
    }

    private static StructuralVariantPosition createSVPosition(long position) {
        return ImmutableStructuralVariantPosition.builder()
                .chromosome(CHROM.chromosome())
                .position(position)
                .id("ID")
                .type(StructuralVariantType.BND)
                .orientation((byte) 1)
                .build();
    }

    private static GenomePosition ratio(final long position) {
        return new GenomePosition() {
            @NotNull
            @Override
            public String chromosome() {
                return CHROM.chromosome();
            }

            @Override
            public long position() {
                return position;
            }
        };
    }

}
