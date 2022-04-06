package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.TestUtils.createDefaultFittedRegion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendLongArmTest
{
    private static final String CHROMOSOME = "13";
    private static final double EPSILON = 1e-10;

    @Test
    public void testCentromereToTelomere()
    {
        final CombinedRegion first = createCombinedRegion(1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion(5001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 5000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        assertCombinedRegion(5001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void testIneligibleChromosome()
    {
        final CombinedRegion first = createCombinedRegion("1", 1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion("1", 5001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 5000, 3, CopyNumberMethod.UNKNOWN, result.get(0));
        assertCombinedRegion(5001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void testExtendThroughStructuralVariantButKeepBreakIntact()
    {

        final CombinedRegion first = createCombinedRegion(1, 1000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion second = createCombinedRegion(1001, 2000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion third = createCombinedRegion(2001, 5000, 3, 0.3, SegmentSupport.BND);
        final CombinedRegion centromere = createCombinedRegion(5001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, third, centromere));
        assertEquals(3, result.size());

        assertCombinedRegion(1, 2000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        assertCombinedRegion(2001, 5000, 2, CopyNumberMethod.LONG_ARM, result.get(1));
        assertCombinedRegion(5001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(2));
    }

    @Test
    public void testDoesNotExtendRight()
    {
        final CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        final CombinedRegion unprocessedRight = createCombinedRegion(20001, 30000, 3, 0.3, SegmentSupport.NONE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(centromere, unprocessedRight));
        assertEquals(2, result.size());

        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(0));
        assertCombinedRegion(20001, 30000, 3, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void extendsToStart()
    {
        final CombinedRegion first = createCombinedRegion(1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion second = createCombinedRegion(5001, 10000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 10000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void extendsThroughUnprocessed()
    {
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
    public void extendsStopsAtProcessed()
    {
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

    @Test
    public void extendsThroughProcessedOnChromosome21()
    {
        final CombinedRegion first = createCombinedRegion("21", 1, 5000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion second = createCombinedRegion("21", 5001, 10000, 3, 0.3, SegmentSupport.NONE);
        final CombinedRegion centromere = createCombinedRegion("21", 10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        first.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, 3);

        final List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 10000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    static void assertCombinedRegion(int start, int end, double expectedCopyNumber, CopyNumberMethod expectedMethod,
            CombinedRegion victim)
    {
        assertEquals(expectedCopyNumber, victim.tumorCopyNumber(), EPSILON);
        assertEquals(expectedMethod, victim.copyNumberMethod());
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
    }

    @NotNull
    private static CombinedRegion createCombinedRegion(String chromosome, int start, int end, double copyNumber, double baf,
            SegmentSupport support)
    {
        final ObservedRegion region = createDefaultFittedRegion(chromosome, start, end);
        region.setTumorCopyNumber(copyNumber);
        region.setTumorBAF(baf);
        region.setSupport(support);
        return new CombinedRegion(region);
    }

    @NotNull
    static CombinedRegion createCombinedRegion(int start, int end, double copyNumber, double baf, SegmentSupport support)
    {
        return createCombinedRegion(CHROMOSOME, start, end, copyNumber, baf, support);
    }
}
