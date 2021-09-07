package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.BND;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.MULTIPLE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.pcf.ImmutablePCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleSegmentFactoryTest
{

    private static final GenomePosition CHROMOSOME_LENGTH = GenomePositions.create("1", 10_000_000);
    private static final GenomePosition CHROMOSOME_CENTROMERE = GenomePositions.create("1", 5_000_001);

    @Test
    public void testEmpty()
    {
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, Collections.emptyList());
        assertEquals(2, segments.size());
        assertPurpleSegment(segments.get(0), 1, CHROMOSOME_CENTROMERE.position() - 1, true, TELOMERE);
        assertPurpleSegment(segments.get(1), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), true, CENTROMERE);
    }

    @Test
    public void testSingleSV()
    {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, 18881).build());
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, clusters);
        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_CENTROMERE.position() - 1, false, BND);
        assertPurpleSegment(segments.get(2), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), false, CENTROMERE);
    }

    @Test
    public void testSVAtCentromere()
    {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, CHROMOSOME_CENTROMERE.position()).build());
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, clusters);
        assertEquals(2, segments.size());

        assertPurpleSegment(segments.get(0), 1, CHROMOSOME_CENTROMERE.position() - 1, true, TELOMERE);
        assertPurpleSegment(segments.get(1), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), false, CENTROMERE);
    }

    @Test
    public void testSingleSVWithRatioSupport()
    {
        final Cluster cluster = addRatios(cluster(17002, 18881), 17050, 19000).build();
        final List<PurpleSegment> segments =
                PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, Lists.newArrayList(cluster));
        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_CENTROMERE.position() - 1, true, BND);
        assertPurpleSegment(segments.get(2), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), true, CENTROMERE);
    }

    @Test
    public void testMultipleSVAtSamePosition()
    {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, 18881).addVariants(variant(18881)).build());
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, clusters);
        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_CENTROMERE.position() - 1, false, MULTIPLE);
        assertPurpleSegment(segments.get(2), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), false, CENTROMERE);
    }

    @Test
    public void testMultipleSVInSameCluster()
    {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, 18881).addVariants(variant(19991)).build());
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, clusters);
        assertEquals(4, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, 19990, false, BND);
        assertPurpleSegment(segments.get(2), 19991, CHROMOSOME_CENTROMERE.position() - 1, false, BND);
        assertPurpleSegment(segments.get(3), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), false, CENTROMERE);
    }

    @Test
    public void testRatiosOnly()
    {
        final Cluster cluster = addRatios(cluster(17002), 18881, 19000).build();
        final List<PurpleSegment> segments =
                PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, Lists.newArrayList(cluster));
        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_CENTROMERE.position() - 1, true, NONE);
        assertPurpleSegment(segments.get(2), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), true, CENTROMERE);
    }

    private static void assertPurpleSegment(@NotNull final PurpleSegment victim, long start, long end, boolean ratioSupport,
            @NotNull final SegmentSupport support)
    {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(ratioSupport, victim.ratioSupport());
        assertEquals(support, victim.support());
    }

    @NotNull
    private static ImmutableCluster.Builder cluster(long start)
    {
        return ImmutableCluster.builder().chromosome(CHROMOSOME_LENGTH.chromosome()).start(start).end(start);
    }

    @NotNull
    private static ImmutableCluster.Builder cluster(long start, long... variants)
    {
        ImmutableCluster.Builder builder = cluster(start);
        for(long position : variants)
        {
            builder.addVariants(variant(position));
        }

        return builder;
    }

    @NotNull
    private static ImmutableCluster.Builder addRatios(@NotNull ImmutableCluster.Builder builder, long... ratios)
    {
        for(long position : ratios)
        {
            builder.addPcfPositions(ImmutablePCFPosition.builder()
                    .chromosome(CHROMOSOME_LENGTH.chromosome())
                    .position(position)
                    .source(PCFSource.TUMOR_RATIO)
                    .minPosition(0)
                    .maxPosition(0)
                    .build());
        }

        return builder;
    }

    @NotNull
    private static SVSegment variant(long position)
    {
        return ImmutableSVSegment.builder()
                .chromosome(CHROMOSOME_LENGTH.chromosome())
                .position(position)
                .type(StructuralVariantType.BND)
                .build();
    }
}
