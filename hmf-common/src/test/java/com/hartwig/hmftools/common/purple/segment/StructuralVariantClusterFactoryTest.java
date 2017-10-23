package com.hartwig.hmftools.common.purple.segment;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.junit.Before;
import org.junit.Test;

public class StructuralVariantClusterFactoryTest {
    private static final String CHROM = "1";

    private StructuralVariantClusterFactory victim;

    @Before
    public void setup() {
        victim = new StructuralVariantClusterFactory(1000);
    }

    @Test
    public void testDefaultClusterBounds() {
        final StructuralVariantPosition sv = createSVPosition(15532);
        final List<StructuralVariantCluster> clusters = victim.cluster(Lists.newArrayList(sv), Collections.emptyList());
        assertEquals(1, clusters.size());
        assertCluster(clusters.get(0), 14001, 17000, 15532);
    }

    @Test
    public void testClusterBoundsWithRatios() {
        final StructuralVariantPosition sv = createSVPosition(15001);
        final List<CobaltRatio> ratios = createRatios();
        final List<StructuralVariantCluster> clusters = victim.cluster(Lists.newArrayList(sv), ratios);
        assertEquals(1, clusters.size());
        assertCluster(clusters.get(0), 12001, 19000, 15001);
    }

    @Test
    public void testTwoSVInsideCluster() {
        final StructuralVariantPosition sv1 = createSVPosition(15532);
        final StructuralVariantPosition sv2 = createSVPosition(17771);
        final List<StructuralVariantCluster> clusters = victim.cluster(Lists.newArrayList(sv1, sv2), Collections.emptyList());
        assertEquals(1, clusters.size());
        assertCluster(clusters.get(0), 14001, 19000, 15532, 17771);
    }

    @Test
    public void testTwoSVOutsideCluster() {
        final StructuralVariantPosition sv1 = createSVPosition(15532);
        final StructuralVariantPosition sv2 = createSVPosition(18881);
        final List<StructuralVariantCluster> clusters = victim.cluster(Lists.newArrayList(sv1, sv2), Collections.emptyList());
        assertEquals(2, clusters.size());
        assertCluster(clusters.get(0), 14001, 17000, 15532);
        assertCluster(clusters.get(1), 17001, 20000, 18881);
    }

    @Test
    public void testTwoSVInsideClusterWithRatio() {
        final StructuralVariantPosition sv1 = createSVPosition(15532);
        final StructuralVariantPosition sv2 = createSVPosition(18881);
        final List<CobaltRatio> ratios = createRatios();
        final List<StructuralVariantCluster> clusters = victim.cluster(Lists.newArrayList(sv1, sv2), ratios);
        assertEquals(1, clusters.size());
        assertCluster(clusters.get(0), 12001, 20000, 15532, 18881);
    }

    private List<CobaltRatio> createRatios() {
        final List<CobaltRatio> ratios = Lists.newArrayList(ratio(11001, true),
                ratio(12001, true),
                ratio(13001, false),
                ratio(14001, false),
                ratio(15001, true),
                ratio(16001, false),
                ratio(17001, false),
                ratio(18001, true));

        return ratios;
    }

    private static void assertCluster(final StructuralVariantCluster cluster, long start, long end, long position) {
        assertEquals(start, cluster.start());
        assertEquals(end, cluster.end());
        assertEquals(position, cluster.firstVariantPosition());
        assertEquals(position, cluster.finalVariantPosition());
        assertEquals(StructuralVariantSupport.BND, cluster.type());
    }

    private static void assertCluster(final StructuralVariantCluster cluster, long start, long end, long firstPosition, long finalPosition) {
        assertEquals(start, cluster.start());
        assertEquals(end, cluster.end());
        assertEquals(firstPosition, cluster.firstVariantPosition());
        assertEquals(finalPosition, cluster.finalVariantPosition());
        assertEquals(StructuralVariantSupport.MULTIPLE, cluster.type());
    }

    private static CobaltRatio ratio(long position, boolean useable) {
        return ratio(position, useable ? 1 : -1);
    }

    private static CobaltRatio ratio(long position, double ratio) {
        return ImmutableCobaltRatio.builder()
                .chromosome(CHROM)
                .position(position)
                .tumorReadCount(0)
                .referenceReadCount(0)
                .referenceGCRatio(1)
                .referenceGCDiploidRatio(1)
                .tumorGCRatio(ratio)
                .build();
    }

    private static StructuralVariantPosition createSVPosition(long position) {
        return ImmutableStructuralVariantPosition.builder()
                .chromosome("CHROM")
                .position(position)
                .id("ID")
                .type(StructuralVariantType.BND)
                .orientation((byte) 1)
                .build();
    }

}
