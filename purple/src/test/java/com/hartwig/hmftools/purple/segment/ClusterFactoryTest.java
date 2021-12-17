package com.hartwig.hmftools.purple.segment;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.utils.pcf.ImmutablePCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Before;
import org.junit.Test;

public class ClusterFactoryTest
{
    private static final String CHROM = "1";

    private ClusterFactory victim;
    private static final int WINDOW = 1000;

    @Before
    public void setup()
    {
        victim = new ClusterFactory(WINDOW);
    }

    @Test
    public void testBoundaries()
    {
        final List<SVSegment> sv = variants(37383599, 37387153);
        final List<PCFPosition> ratios = createRatioBreaks(36965001, 37381001, 37382001, 37384001, 37387001, 37389001);
        final List<CobaltRatio> cobalt = cobalt(37380001, true, false, true, true, true, true, true);

        final List<Cluster> clusters = victim.cluster(sv, ratios, cobalt);
        assertEquals(4, clusters.size());
        assertRatioInCluster(clusters.get(0), 36964002, 36965001);
        assertRatioInCluster(clusters.get(1), 37380002, 37381001, 37382001);
        assertCluster(clusters.get(2), 37382002, 37383599, 37383599, 37384001, 37384001);
        assertCluster(clusters.get(3), 37386002, 37387153, 37387153, 37387001, 37389001);
    }

    @Test
    public void testWindowStartWithRatios()
    {
        final List<CobaltRatio> cobalt = cobalt(37380001, true, false, true, true, true, true, true);

        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37381001, 6, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37381002, 6, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37381999, 6, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37382000, 6, cobalt));

        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37382001, 6, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37382002, 6, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37382999, 6, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37383000, 6, cobalt));

        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37383001, 6, cobalt));
        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37383002, 6, cobalt));
        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37383999, 6, cobalt));
        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37384000, 6, cobalt));
    }

    @Test
    public void testWindowStartWithoutRatios()
    {
        final List<CobaltRatio> cobalt = Lists.newArrayList();

        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37381001, -1, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37381002, -1, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37381999, -1, cobalt));
        assertEquals(37380002, victim.earliestDetectableCopyNumberChangePosition(37382000, -1, cobalt));

        assertEquals(37381002, victim.earliestDetectableCopyNumberChangePosition(37382001, -1, cobalt));
        assertEquals(37381002, victim.earliestDetectableCopyNumberChangePosition(37382002, -1, cobalt));
        assertEquals(37381002, victim.earliestDetectableCopyNumberChangePosition(37382999, -1, cobalt));
        assertEquals(37381002, victim.earliestDetectableCopyNumberChangePosition(37383000, -1, cobalt));

        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37383001, -1, cobalt));
        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37383002, -1, cobalt));
        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37383999, -1, cobalt));
        assertEquals(37382002, victim.earliestDetectableCopyNumberChangePosition(37384000, -1, cobalt));
    }

    @Test
    public void testDefaultClusterBounds()
    {
        final SVSegment sv = createSVPosition(15532);
        final List<Cluster> clusters = victim.cluster(Lists.newArrayList(sv), Collections.emptyList(), Collections.emptyList());
        assertEquals(1, clusters.size());
        assertVariantInCluster(clusters.get(0), 14002, 15532);
    }

    @Test
    public void testClusterBoundsWithRatios()
    {
        final List<SVSegment> sv = variants(15532);
        final List<CobaltRatio> ratios = createRatios();
        final List<Cluster> clusters = victim.cluster(sv, Collections.emptyList(), ratios);
        assertEquals(1, clusters.size());
        assertVariantInCluster(clusters.get(0), 12002, 15532);
    }

    @Test
    public void testTwoSVInsideCluster()
    {
        final List<SVSegment> sv = variants(15532, 16771);
        final List<Cluster> clusters = victim.cluster(sv, Collections.emptyList(), Collections.emptyList());
        assertEquals(1, clusters.size());
        assertVariantsInCluster(clusters.get(0), 14002, 15532, 16771);
    }

    @Test
    public void testTwoSVOutsideCluster()
    {
        final List<SVSegment> sv = variants(15532, 17881);
        final List<Cluster> clusters = victim.cluster(sv, Collections.emptyList(), Collections.emptyList());
        assertEquals(2, clusters.size());
        assertVariantInCluster(clusters.get(0), 14002, 15532);
        assertVariantInCluster(clusters.get(1), 16002, 17881);
    }

    @Test
    public void testTwoSVInsideClusterWithRatio()
    {
        final List<SVSegment> sv = variants(15532, 18881);
        final List<CobaltRatio> ratios = createRatios();
        final List<Cluster> clusters = victim.cluster(sv, Collections.emptyList(), ratios);
        assertEquals(1, clusters.size());
        assertVariantsInCluster(clusters.get(0), 12002, 15532, 18881);
    }

    private static void assertRatioInCluster(final Cluster cluster, int start, int position)
    {
        assertCluster(cluster, start, null, null, position, position);
    }

    private static void assertRatioInCluster(final Cluster cluster, int start, int firstPosition, int finalPosition)
    {
        assertCluster(cluster, start, null, null, firstPosition, finalPosition);
    }

    private static void assertVariantInCluster(final Cluster cluster, int start, int position)
    {
        assertCluster(cluster, start, position, position, null, null);
    }

    private static void assertVariantsInCluster(final Cluster cluster, int start, int firstPosition, int finalPosition)
    {
        assertCluster(cluster, start, firstPosition, finalPosition, null, null);
    }

    private static void assertCluster(@NotNull final Cluster cluster, int start, @Nullable Integer firstVariant, @Nullable Integer finalVariant,
            @Nullable Integer firstRatio, @Nullable Integer finalRatio)
    {
        assertEquals(start, cluster.start());
        assertEquals(Math.max(finalVariant == null ? 0 : finalVariant, finalRatio == null ? 0 : finalRatio), cluster.end());
        assertEquals(firstVariant, firstVariant(cluster));
        assertEquals(finalVariant, finalVariant(cluster));
        assertEquals(firstRatio, firstRatio(cluster));
        assertEquals(finalRatio, finalRatio(cluster));
    }

    @NotNull
    private static List<CobaltRatio> createRatios()
    {
        return cobalt(11001, true, true, false, false, true, false, false, true);
    }

    @NotNull
    private static CobaltRatio cobalt(int position, boolean useable)
    {
        return ratio(position, useable ? 1 : -1);
    }

    @NotNull
    private static CobaltRatio ratio(int position, double ratio)
    {
        return PurpleTestUtils.cobalt(CHROM, position, ratio).build();
    }

    @NotNull
    private static SVSegment createSVPosition(int position)
    {
        return ImmutableSVSegment.builder().chromosome(CHROM).position(position).type(StructuralVariantType.BND).build();
    }

    @NotNull
    private static List<CobaltRatio> cobalt(int startPosition, boolean... usable)
    {
        final List<CobaltRatio> result = Lists.newArrayList();
        int offset = 0;
        for(boolean isUsable : usable)
        {
            result.add(cobalt(startPosition + offset, isUsable));
            offset += WINDOW;
        }
        return result;
    }

    @NotNull
    private static List<SVSegment> variants(int... positions)
    {
        final List<SVSegment> result = Lists.newArrayList();
        for(int position : positions)
        {
            result.add(createSVPosition(position));
        }

        return result;
    }

    @NotNull
    private static List<PCFPosition> createRatioBreaks(int... positions)
    {
        final List<PCFPosition> result = Lists.newArrayList();
        for(int position : positions)
        {
            result.add(ratio(position));
        }

        return result;
    }

    @NotNull
    private static PCFPosition ratio(int position)
    {
        return ImmutablePCFPosition.builder()
                .chromosome(CHROM)
                .position(position)
                .source(PCFSource.TUMOR_RATIO)
                .minPosition(0)
                .maxPosition(0)
                .build();
    }

    @Nullable
    private static Integer firstVariant(@NotNull Cluster cluster)
    {
        return cluster.variants().isEmpty() ? null : cluster.variants().get(0).position();
    }

    @Nullable
    private static Integer finalVariant(@NotNull Cluster cluster)
    {
        return cluster.variants().isEmpty() ? null : cluster.variants().get(cluster.variants().size() - 1).position();
    }

    @Nullable
    private static Integer firstRatio(@NotNull Cluster cluster)
    {
        return cluster.ratios().isEmpty() ? null : cluster.ratios().get(0).position();
    }

    @Nullable
    private static Integer finalRatio(@NotNull Cluster cluster)
    {
        return cluster.ratios().isEmpty() ? null : cluster.ratios().get(cluster.ratios().size() - 1).position();
    }
}
