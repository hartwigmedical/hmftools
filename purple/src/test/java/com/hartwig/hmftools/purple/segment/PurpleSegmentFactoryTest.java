package com.hartwig.hmftools.purple.segment;

import static com.hartwig.hmftools.common.purple.SegmentSupport.BND;
import static com.hartwig.hmftools.common.purple.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.MULTIPLE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.TELOMERE;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
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
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, 18881));
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, clusters);
        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_CENTROMERE.position() - 1, false, BND);
        assertPurpleSegment(segments.get(2), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), false, CENTROMERE);
    }

    @Test
    public void testSVAtCentromere()
    {
        final List<Cluster> clusters = Lists.newArrayList(cluster(17001, CHROMOSOME_CENTROMERE.position()));
        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, clusters);
        assertEquals(2, segments.size());

        assertPurpleSegment(segments.get(0), 1, CHROMOSOME_CENTROMERE.position() - 1, true, TELOMERE);
        assertPurpleSegment(segments.get(1), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), false, CENTROMERE);
    }

    @Test
    public void testSingleSVWithRatioSupport()
    {
        final Cluster cluster = cluster(17002, 18881);
        addRatios(cluster, 17050, 19000);

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
        Cluster cluster = cluster(17001, 18881);
        cluster.Variants.add(variant(18881));

        final List<Cluster> clusters = Lists.newArrayList(cluster);

        final List<PurpleSegment> segments = PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, clusters);

        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_CENTROMERE.position() - 1, false, MULTIPLE);
        assertPurpleSegment(segments.get(2), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), false, CENTROMERE);
    }

    @Test
    public void testMultipleSVInSameCluster()
    {
        Cluster cluster = cluster(17001, 18881);
        cluster.Variants.add(variant(19991));
        final List<Cluster> clusters = Lists.newArrayList(cluster);

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
        final Cluster cluster = cluster(17002);
        addRatios(cluster, 18881, 19000);

        final List<PurpleSegment> segments =
                PurpleSegmentFactory.create(CHROMOSOME_CENTROMERE, CHROMOSOME_LENGTH, Lists.newArrayList(cluster));

        assertEquals(3, segments.size());
        assertPurpleSegment(segments.get(0), 1, 18880, true, TELOMERE);
        assertPurpleSegment(segments.get(1), 18881, CHROMOSOME_CENTROMERE.position() - 1, true, NONE);
        assertPurpleSegment(segments.get(2), CHROMOSOME_CENTROMERE.position(), CHROMOSOME_LENGTH.position(), true, CENTROMERE);
    }

    private static void assertPurpleSegment(final PurpleSegment victim, int start, int end, boolean ratioSupport, final SegmentSupport support)
    {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(ratioSupport, victim.RatioSupport);
        assertEquals(support, victim.Support);
    }

    @NotNull
    private static Cluster cluster(int start)
    {
        return new Cluster(CHROMOSOME_LENGTH.chromosome(), start, start);
    }

    @NotNull
    private static Cluster cluster(int start, int... variants)
    {
        Cluster cluster = cluster(start);
        for(int position : variants)
        {
            cluster.Variants.add(variant(position));
        }

        return cluster;
    }

    @NotNull
    private static void addRatios(final Cluster cluster, int... ratios)
    {
        for(int position : ratios)
        {
            cluster.PcfPositions.add(new PCFPosition(PCFSource.TUMOR_RATIO, CHROMOSOME_LENGTH.chromosome(), position));
        }
    }

    @NotNull
    private static SVSegment variant(int position)
    {
        return new SVSegment(CHROMOSOME_LENGTH.chromosome(), position, StructuralVariantType.BND);
    }
}
