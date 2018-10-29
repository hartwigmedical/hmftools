package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

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


    private void assertInferRegion(@NotNull final ExtendDiploidBAF.InferRegion victim, int expectedLeftSource, int expectedLeftTarget,
            int expectedRightTarget, int expectedRightSource) {
        assertEquals(expectedLeftSource, victim.leftSource);
        assertEquals(expectedLeftTarget, victim.leftTarget);
        assertEquals(expectedRightTarget, victim.rightTarget);
        assertEquals(expectedRightSource, victim.rightSource);
    }

    @NotNull
    private static CombinedRegion create(long start, long end, SegmentSupport support, int bafCount) {
        return new BafWeightedRegion(PurpleDatamodelTest.createDefaultFittedRegion("1", start, end)
                .support(support)
                .bafCount(bafCount)
                .build());
    }

    @NotNull
    private static ExtendDiploidBAF.InferRegion inferRegion(final int leftSource, final int leftTarget, final int rightTarget, final int rightSource) {
        return new ExtendDiploidBAF.InferRegion(leftSource, leftTarget, rightTarget, rightSource);
    }

}
