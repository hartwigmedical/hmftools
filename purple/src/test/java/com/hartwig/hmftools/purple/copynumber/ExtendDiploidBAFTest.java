package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.MiscTestUtils.createDefaultFittedRegion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendDiploidBAFTest
{
    private static final double EPSILON = 1e-10;

    @Test
    public void testDUPWithLOH()
    {
        final Map<Integer,Integer> dupMap = Maps.newHashMap();
        dupMap.put(1001, 2001);
        final ExtendDiploidBAF victim = new ExtendDiploidBAF(dupMap);

        CombinedRegion leftWithLOH = create(1, 1000, SegmentSupport.TELOMERE, 10, 1, 2);
        CombinedRegion dup = create(1001, 2000, SegmentSupport.DUP, 0, 1, 4);
        CombinedRegion del = create(1001, 2000, SegmentSupport.DEL, 0, 1, 4);
        CombinedRegion rightWithLOH = create(2001, 3000, SegmentSupport.DUP, 1, 1, 2);
        CombinedRegion rightWithoutLOH = create(2001, 3000, SegmentSupport.DUP, 1, 0.5, 2);
        ExtendDiploidBAF.InferRegion inferRegion = new ExtendDiploidBAF.InferRegion(0, 1, 1, 2);

        assertTrue(victim.isSimpleDupSurroundedByLOH(inferRegion, Lists.newArrayList(leftWithLOH, dup, rightWithLOH)));
        assertFalse(victim.isSimpleDupSurroundedByLOH(inferRegion, Lists.newArrayList(leftWithLOH, dup, rightWithoutLOH)));
        assertFalse(victim.isSimpleDupSurroundedByLOH(inferRegion, Lists.newArrayList(leftWithLOH, del, rightWithLOH)));
    }

    @Test
    public void testTargetAllele()
    {
        assertAlleleTarget(1, 0.5, 0.5, 0);
        assertAlleleTarget(1, 1, 1, 0);
        assertAlleleTarget(1, 2, 1, 1);
        assertAlleleTarget(1, 2.5, 1.5, 1);
        assertAlleleTarget(1, 3, 2, 1);
    }

    private void assertAlleleTarget(double majorTarget, double copyNumber, double expectedMajor, double expectedMinor)
    {
        double baf = ExtendDiploidBAF.bafForTargetAllele(majorTarget, copyNumber);
        double actualMajor = baf * copyNumber;
        double actualMinor = (1 - baf) * copyNumber;
        assertEquals(expectedMajor, actualMajor, EPSILON);
        assertEquals(expectedMinor, actualMinor, EPSILON);
    }

    @Test
    public void testShortArm()
    {
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
    public void testTinyArmWithMultipleSources()
    {
        CombinedRegion cr1 = create(1, 1000, SegmentSupport.NONE, 1000, 1, 1.9);
        CombinedRegion cr2 = create(1001, 1030, SegmentSupport.NONE, 0, 0, 2.2);
        CombinedRegion cr2a = create(1001, 1031, SegmentSupport.NONE, 0, 0, 2.2);
        CombinedRegion cr3 = create(1031, 2000, SegmentSupport.NONE, 100, 1, 1.9);
        CombinedRegion cr3a = create(1032, 2000, SegmentSupport.NONE, 100, 1, 1.9);
        CombinedRegion cr4 = create(2001, 3000, SegmentSupport.NONE, 100, 0.666, 3);

        List<CombinedRegion> result = new ExtendDiploidBAF(Collections.emptyList()).extendBAF(Lists.newArrayList(cr1, cr2, cr3, cr4));
        assertEquals(0, result.get(1).region().minorAlleleCopyNumber(), EPSILON);

        result = new ExtendDiploidBAF(Collections.emptyList()).extendBAF(Lists.newArrayList(cr1, cr2a, cr3a, cr4));
        assertEquals(0.3, result.get(1).region().minorAlleleCopyNumber(), EPSILON);
    }

    @Test
    public void testSmallRegionWithinLargeLOHRegion()
    {
        CombinedRegion cr1 = create(80_000_000, 100_000_000, SegmentSupport.NONE, 1000, 0.66, 3);
        CombinedRegion cr2 = create(100_000_001, 101_000_000, SegmentSupport.NONE, 1000, 1, 2.0);
        CombinedRegion cr3 = create(101_000_001, 101_001_000, SegmentSupport.NONE, 0, 0, 3.0);
        CombinedRegion cr4 = create(101_001_001, 104_000_000, SegmentSupport.NONE, 1000, 1, 2.0);

        List<CombinedRegion> result = new ExtendDiploidBAF(Collections.emptyList()).extendBAF(Lists.newArrayList(cr1, cr2, cr3, cr4));
        assertEquals(0, result.get(2).region().minorAlleleCopyNumber(), EPSILON);
    }

    @Test
    public void testSmallRegionWithinLargeLOHRegionAfter()
    {
        CombinedRegion cr1 = create(100_000_001, 101_000_000, SegmentSupport.NONE, 1000, 1, 2.0);
        CombinedRegion cr2 = create(101_000_001, 101_001_000, SegmentSupport.NONE, 0, 0, 3.0);
        CombinedRegion cr3 = create(101_001_001, 102_001_000, SegmentSupport.NONE, 1000, 1, 2.0);
        CombinedRegion cr4 = create(102_001_001, 104_000_000, SegmentSupport.NONE, 1000, 0.66, 3);

        List<CombinedRegion> result = new ExtendDiploidBAF(Collections.emptyList()).extendBAF(Lists.newArrayList(cr1, cr2, cr3, cr4));
        assertEquals(0, result.get(1).region().minorAlleleCopyNumber(), EPSILON);
    }

    @Test
    public void testMinorOrMajorMovedTargetPloidyWithCommonMinor()
    {
        final ObservedRegion left = create(1.01, 2.01);
        final ObservedRegion right = create(1.02, 3.02);
        final ObservedRegion sourceLikeLeft = create(1.03, 2.03);
        final ObservedRegion sourceLikeRight = create(1.04, 3.04);

        assertEquals(1.01, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(left, left, right), EPSILON);
        assertEquals(1.02, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(right, left, right), EPSILON);

        assertEquals(1.03, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeLeft, left, right), EPSILON);
        assertEquals(1.04, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeRight, left, right), EPSILON);
    }

    @Test
    public void testMinorOrMajorMovedTargetPloidyWithCommonMajor()
    {
        final ObservedRegion left = create(1.01, 2.01);
        final ObservedRegion right = create(0.02, 2.02);
        final ObservedRegion sourceLikeLeft = create(1.03, 2.03);
        final ObservedRegion sourceLikeRight = create(0.04, 2.04);

        assertEquals(2.01, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(left, left, right), EPSILON);
        assertEquals(2.02, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(right, left, right), EPSILON);

        assertEquals(2.03, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeLeft, left, right), EPSILON);
        assertEquals(2.04, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeRight, left, right), EPSILON);
    }

    @Test
    public void testMinorOrMajorMovedTargetPloidyWithMajorMinorCross()
    {
        final ObservedRegion left = create(1.01, 2.01);
        final ObservedRegion right = create(2.02, 3.02);
        final ObservedRegion sourceLikeLeft = create(1.03, 2.03);
        final ObservedRegion sourceLikeRight = create(2.04, 3.04);

        assertEquals(2.01, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(left, left, right), EPSILON);
        assertEquals(2.02, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(right, left, right), EPSILON);

        assertEquals(2.03, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeLeft, left, right), EPSILON);
        assertEquals(2.04, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeRight, left, right), EPSILON);
    }

    @Test
    public void testMinorOrMajorMovedTargetPloidyWithMinorMajorCross()
    {
        final ObservedRegion left = create(1.01, 2.01);
        final ObservedRegion right = create(0.02, 1.02);
        final ObservedRegion sourceLikeLeft = create(0.03, 1.03);
        final ObservedRegion sourceLikeRight = create(0.04, 1.04);

        assertEquals(1.01, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(left, left, right), EPSILON);
        assertEquals(1.02, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(right, left, right), EPSILON);

        assertEquals(1.03, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeLeft, left, right), EPSILON);
        assertEquals(1.04, ExtendDiploidBAF.minorOrMajorMovedTargetPloidy(sourceLikeRight, left, right), EPSILON);
    }

    private void assertInferRegion(@NotNull final ExtendDiploidBAF.InferRegion victim, int expectedLeftSource, int expectedLeftTarget,
            int expectedRightTarget, int expectedRightSource)
    {
        assertEquals(expectedLeftSource, victim.leftSourceIndex);
        assertEquals(expectedLeftTarget, victim.leftTargetIndex);
        assertEquals(expectedRightTarget, victim.rightTargetIndex);
        assertEquals(expectedRightSource, victim.rightSourceIndex);
    }

    @NotNull
    private static CombinedRegion create(int start, int end, SegmentSupport support, int bafCount)
    {
        ObservedRegion region = createDefaultFittedRegion("1", start, end);
        region.setSupport(support);
        region.setBafCount(bafCount);
        return new CombinedRegion(region);
    }

    @NotNull
    private static CombinedRegion create(int start, int end, SegmentSupport support, int bafCount, double baf, double copyNumber)
    {
        ObservedRegion region = createDefaultFittedRegion("1", start, end);
        region.setSupport(support);
        region.setBafCount(bafCount);
        region.setTumorBAF(baf);
        region.setTumorCopyNumber(copyNumber);
        return new CombinedRegion(region);
    }

    @NotNull
    private static ObservedRegion create(double minorAllele, double majorAllele)
    {
        double copyNumber = minorAllele + majorAllele;
        double baf = majorAllele / copyNumber;
        ObservedRegion region = createDefaultFittedRegion("1", 1, 1000);
        region.setTumorBAF(baf);
        region.setTumorCopyNumber(copyNumber);
        return region;
    }
}
