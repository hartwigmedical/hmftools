package com.hartwig.hmftools.common.purple.segment;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.BND;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.MULTIPLE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.chromosome.ImmutableChromosomeLength;
import com.hartwig.hmftools.common.pcf.ImmutablePCFPosition;
import com.hartwig.hmftools.common.pcf.PCFSource;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleSegmentFactoryTest {
    static final ChromosomeLength CHROMOSOME_LENGTH = ImmutableChromosomeLength.builder().chromosome("1").position(10_000_000).build();

    @Test
    public void testEmpty() {
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_LENGTH, Collections.emptyList());
        assertEquals(1, segments.size());
        assertPurpleSegment(segments.get(0), 1, CHROMOSOME_LENGTH.position(), true, TELOMERE);
    }

    @Test
    public void testSingleSV() {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, 18881).build());
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_LENGTH, clusters);
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_LENGTH.position(), false, BND);
    }

    @Test
    public void testSingleSVWithRatioSupport() {
        final Cluster cluster = addRatios(cluster(17002, 18881), 17050, 19000).build();
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_LENGTH, Lists.newArrayList(cluster));
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_LENGTH.position(), true, BND);
    }


    @Test
    public void testMultipleSVAtSamePosition() {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, 18881).addVariants(variant(18881)).build());
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_LENGTH, clusters);
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_LENGTH.position(), false, MULTIPLE);
    }

    @Test
    public void testMultipleSVInSameCluster() {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, 18881).addVariants(variant(19991)).build());
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_LENGTH, clusters);
        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, 19990, false, BND);
        assertPurpleSegment(segments.get(2), 19991, CHROMOSOME_LENGTH.position(), false, BND);
    }


    @Test
    public void testRatiosOnly() {
        final Cluster cluster = addRatios(cluster(17002), 18881, 19000).build();
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_LENGTH, Lists.newArrayList(cluster));
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_LENGTH.position(), true, NONE);
    }

    private static void assertPurpleSegment(@NotNull final PurpleSegment victim, long start, long end, boolean ratioSupport,
            @NotNull final SegmentSupport support) {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(ratioSupport, victim.ratioSupport());
        assertEquals(support, victim.support());
    }

    private static ImmutableCluster.Builder cluster(long start) {
        return ImmutableCluster.builder().chromosome(CHROMOSOME_LENGTH.chromosome()).start(start).end(start);
    }

    private static ImmutableCluster.Builder cluster(long start, long... variants) {

        ImmutableCluster.Builder builder = cluster(start);
        for (long position : variants) {
            builder.addVariants(variant(position));
        }

        return builder;
    }

    private static ImmutableCluster.Builder addRatios(ImmutableCluster.Builder builder, long... ratios) {
        for (long position : ratios) {
            builder.addPcfPositions(ImmutablePCFPosition.builder()
                    .chromosome(CHROMOSOME_LENGTH.chromosome())
                    .position(position)
                    .source(PCFSource.TUMOR_RATIO)
                    .build());
        }

        return builder;
    }

    private static StructuralVariantPosition variant(long position) {
        return ImmutableStructuralVariantPosition.builder()
                .chromosome(CHROMOSOME_LENGTH.chromosome())
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
                return CHROMOSOME_LENGTH.chromosome();
            }

            @Override
            public long position() {
                return position;
            }
        };
    }

}
