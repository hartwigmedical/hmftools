package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendLongArmTest {

    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1e-10;

    @Test
    public void testCentromereToTelomere() {
        final CombinedRegion first = createCombinedRegion(1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion(5001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 5000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        assertCombinedRegion(5001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void testDoesNotExtendRight() {
        final CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        final CombinedRegion unprocessedRight = createCombinedRegion(20001, 30000, 3, 0.3, SegmentSupport.NONE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(centromere, unprocessedRight));
        assertEquals(2, result.size());

        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(0));
        assertCombinedRegion(20001, 30000, 3, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void extendsToStart() {
        final CombinedRegion first = createCombinedRegion(1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion second = createCombinedRegion(5001, 10000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 10000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void extendsThroughUnprocessed() {
        final CombinedRegion first = createCombinedRegion(1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion second = createCombinedRegion(5001, 10000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        first.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, 3);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(3, result.size());

        assertCombinedRegion(1, 5000, 3, CopyNumberMethod.STRUCTURAL_VARIANT, result.get(0));
        assertCombinedRegion(5001, 10000, 2, CopyNumberMethod.LONG_ARM, result.get(1));
        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(2));
    }

    @Test
    public void extendsStopsAtProcessed() {
        final CombinedRegion first = createCombinedRegion(1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion second = createCombinedRegion(5001, 10000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        second.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, 3);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(3, result.size());

        assertCombinedRegion(1, 5000, 3, CopyNumberMethod.UNKNOWN, result.get(0));
        assertCombinedRegion(5001, 10000, 3, CopyNumberMethod.STRUCTURAL_VARIANT, result.get(1));
        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(2));
    }

    static void assertCombinedRegion(long start, long end, double expectedCopyNumber, CopyNumberMethod expectedMethod,
            CombinedRegion victim) {
        assertEquals(expectedCopyNumber, victim.tumorCopyNumber(), EPSILON);
        assertEquals(expectedMethod, victim.copyNumberMethod());
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
    }

    @NotNull
    static CombinedRegion createCombinedRegion(long start, long end, double copyNumber, double baf, SegmentSupport support) {
        final FittedRegion region = PurpleDatamodelTest.createDefaultFittedRegion(CHROMOSOME, start, end)
                .tumorCopyNumber(copyNumber)
                .tumorBAF(baf)
                .support(support)
                .build();

        return new CombinedRegion(region);
    }
}
