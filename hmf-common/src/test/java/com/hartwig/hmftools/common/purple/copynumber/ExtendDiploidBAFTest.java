package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendDiploidBAFTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testTargetAllele() {
        assertAlleleTarget(1, 0.5, 0.5, 0);
        assertAlleleTarget(1, 1, 1, 0);
        assertAlleleTarget(1, 2, 1, 1);
        assertAlleleTarget(1, 2.5, 1.5, 1);
        assertAlleleTarget(1, 3, 2, 1);
    }

    private void assertAlleleTarget(double majorTarget, double copyNumber,  double expectedMajor, double expectedMinor) {
        double baf = ExtendDiploidBAF.bafForTargetAllele(majorTarget, copyNumber);
        double actualMajor = baf * copyNumber;
        double actualMinor = (1-baf) * copyNumber;
        assertEquals(expectedMajor, actualMajor, EPSILON);
        assertEquals(expectedMinor, actualMinor, EPSILON);
    }

    @Test
    public void testShortArm() {
        CombinedRegion cr1 = create(1, 1000, SegmentSupport.NONE, 0);
        CombinedRegion cr2 = create(1001, 2000, SegmentSupport.NONE, 0);
        CombinedRegion cr3 = create(2001, 3000, SegmentSupport.CENTROMERE, 100);
        CombinedRegion cr4 = create(3001, 4000, SegmentSupport.NONE, 0);
        CombinedRegion cr5 = create(4001, 5000, SegmentSupport.NONE, 1);

        assertInferRegion(ExtendDiploidBAF.nextRegion(true, Lists.newArrayList(cr1, cr2, cr3, cr4, cr5)), -1, 0, 1, 2);
        assertInferRegion(ExtendDiploidBAF.nextRegion(false, Lists.newArrayList(cr1, cr2, cr3, cr4, cr5)), 2, 3, 3, 4);
        assertInferRegion(ExtendDiploidBAF.nextRegion(false, Lists.newArrayList(cr1, cr2, cr3, cr4)), 2, 3, 3, -1);

    }

    @Test
    public void testTinyArmWithMultipleSources() {
        CombinedRegion cr1 = create(1, 1000, SegmentSupport.NONE, 1000, 1, 1.9);
        CombinedRegion cr2 = create(1001, 1030, SegmentSupport.NONE, 0, 0, 2.2);
        CombinedRegion cr2a = create(1001, 1031, SegmentSupport.NONE, 0, 0, 2.2);
        CombinedRegion cr3 = create(1031, 2000, SegmentSupport.NONE, 100, 1, 1.9);
        CombinedRegion cr3a = create(1032, 2000, SegmentSupport.NONE, 100, 1, 1.9);
        CombinedRegion cr4 = create(2001, 3000, SegmentSupport.NONE, 100, 0.666, 3);


        List<CombinedRegion> result =  ExtendDiploidBAF.extendBAF(Lists.newArrayList(cr1, cr2, cr3, cr4));
        assertEquals(0, result.get(1).region().minorAllelePloidy(),  EPSILON);

        result =  ExtendDiploidBAF.extendBAF(Lists.newArrayList(cr1, cr2a, cr3a, cr4));
        assertEquals(0.3, result.get(1).region().minorAllelePloidy(),  EPSILON);
    }


    @Test
    public void testSmallRegionWithinLargeLOHRegion() {
        CombinedRegion cr1 = create(80_000_000, 100_000_000, SegmentSupport.NONE, 1000, 0.66, 3);
        CombinedRegion cr2 = create(100_000_001, 101_000_000, SegmentSupport.NONE, 1000, 1, 2.0);
        CombinedRegion cr3 = create(101_000_001, 101_001_000, SegmentSupport.NONE, 0, 0, 3.0);
        CombinedRegion cr4 = create(101_001_001, 104_000_000, SegmentSupport.NONE, 1000, 1, 2.0);

        List<CombinedRegion> result =  ExtendDiploidBAF.extendBAF(Lists.newArrayList(cr1, cr2, cr3, cr4));
        assertEquals(0, result.get(2).region().minorAllelePloidy(),  EPSILON);
    }

    @Test
    public void testSmallRegionWithinLargeLOHRegionAfter() {
        CombinedRegion cr1 = create(100_000_001, 101_000_000, SegmentSupport.NONE, 1000, 1, 2.0);
        CombinedRegion cr2 = create(101_000_001, 101_001_000, SegmentSupport.NONE, 0, 0, 3.0);
        CombinedRegion cr3 = create(101_001_001, 102_001_000, SegmentSupport.NONE, 1000, 1, 2.0);
        CombinedRegion cr4 = create(102_001_001, 104_000_000, SegmentSupport.NONE, 1000, 0.66, 3);

        List<CombinedRegion> result =  ExtendDiploidBAF.extendBAF(Lists.newArrayList(cr1, cr2, cr3, cr4));
        assertEquals(0, result.get(1).region().minorAllelePloidy(),  EPSILON);
    }

    private void assertInferRegion(@NotNull final ExtendDiploidBAF.InferRegion victim, int expectedLeftSource, int expectedLeftTarget,
            int expectedRightTarget, int expectedRightSource) {
        assertEquals(expectedLeftSource, victim.leftSourceIndex);
        assertEquals(expectedLeftTarget, victim.leftTargetIndex);
        assertEquals(expectedRightTarget, victim.rightTargetIndex);
        assertEquals(expectedRightSource, victim.rightSourceIndex);
    }

    @NotNull
    private static CombinedRegion create(long start, long end, SegmentSupport support, int bafCount) {
        return new BafWeightedRegion(PurpleDatamodelTest.createDefaultFittedRegion("1", start, end)
                .support(support)
                .bafCount(bafCount)
                .build());
    }

    @NotNull
    private static CombinedRegion create(long start, long end, SegmentSupport support, int bafCount, double baf, double copyNumber) {
        return new BafWeightedRegion(PurpleDatamodelTest.createDefaultFittedRegion("1", start, end)
                .support(support)
                .bafCount(bafCount)
                .tumorBAF(baf)
                .tumorCopyNumber(copyNumber)
                .build());
    }
}
