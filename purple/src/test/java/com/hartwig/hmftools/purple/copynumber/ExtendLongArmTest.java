package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.MiscTestUtils.createDefaultFittedRegion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.junit.Test;

public class ExtendLongArmTest
{
    private static final String CHROMOSOME_NORMAL = "12";
    private static final String CHROMOSOME_SHORT = "13";
    private static final double EPSILON = 1e-10;

    @Test
    public void testCentromereToTelomere()
    {
        CombinedRegion first = createCombinedRegionShortArm(1, 5000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion centromere = createCombinedRegionShortArm(5001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        centromere.setTumorCopyNumber(CopyNumberMethod.BAF_WEIGHTED, 2);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 5000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        // assertCombinedRegion(5001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void testIneligibleChromosome()
    {
        CombinedRegion first = createCombinedRegion("1", 1, 5000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion centromere = createCombinedRegion("1", 5001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 5000, 3, CopyNumberMethod.UNKNOWN, result.get(0));
        assertCombinedRegion(5001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void testExtendThroughStructuralVariantButKeepBreakIntact()
    {
        CombinedRegion first = createCombinedRegionShortArm(1, 1000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion second = createCombinedRegionShortArm(1001, 2000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion third = createCombinedRegionShortArm(2001, 5000, 3, 0.3, SegmentSupport.BND);
        CombinedRegion centromere = createCombinedRegionShortArm(5001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        centromere.setTumorCopyNumber(CopyNumberMethod.BAF_WEIGHTED, 2);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, third, centromere));
        assertEquals(3, result.size());

        assertCombinedRegion(1, 2000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        assertCombinedRegion(2001, 5000, 2, CopyNumberMethod.LONG_ARM, result.get(1));
        // assertCombinedRegion(5001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(2));
    }

    @Test
    public void testDoesNotExtendRight()
    {
        CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        CombinedRegion unprocessedRight = createCombinedRegion(20001, 30000, 3, 0.3, SegmentSupport.NONE);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(centromere, unprocessedRight));
        assertEquals(2, result.size());

        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(0));
        assertCombinedRegion(20001, 30000, 3, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void extendsToStart()
    {
        CombinedRegion first = createCombinedRegionShortArm(1, 5000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion second = createCombinedRegionShortArm(5001, 10000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion centromere = createCombinedRegionShortArm(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        centromere.setTumorCopyNumber(CopyNumberMethod.BAF_WEIGHTED, 2);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 10000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        // assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    @Test
    public void extendsThroughUnprocessed()
    {
        CombinedRegion first = createCombinedRegion("Y", 1, 5000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion second = createCombinedRegion("Y",5001, 10000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion centromere = createCombinedRegion("Y",10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        first.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, 3);
        centromere.setTumorCopyNumber(CopyNumberMethod.BAF_WEIGHTED, 2);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(3, result.size());

        assertCombinedRegion(1, 5000, 3, CopyNumberMethod.STRUCTURAL_VARIANT, result.get(0));
        assertCombinedRegion(5001, 10000, 2, CopyNumberMethod.LONG_ARM, result.get(1));
        // assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(2));
    }

    @Test
    public void extendsStopsAtProcessed()
    {
        CombinedRegion first = createCombinedRegion(1, 5000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion second = createCombinedRegion(5001, 10000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion centromere = createCombinedRegion(10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        second.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, 3);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(3, result.size());

        assertCombinedRegion(1, 5000, 3, CopyNumberMethod.UNKNOWN, result.get(0));
        assertCombinedRegion(5001, 10000, 3, CopyNumberMethod.STRUCTURAL_VARIANT, result.get(1));
        assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(2));
    }

    @Test
    public void extendsThroughProcessedOnChromosome21()
    {
        CombinedRegion first = createCombinedRegion("21", 1, 5000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion second = createCombinedRegion("21", 5001, 10000, 3, 0.3, SegmentSupport.NONE);
        CombinedRegion centromere = createCombinedRegion("21", 10001, 20000, 2, 0.5, SegmentSupport.CENTROMERE);
        first.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, 3);
        centromere.setTumorCopyNumber(CopyNumberMethod.BAF_WEIGHTED, 2);

        List<CombinedRegion> result = ExtendLongArm.extendLongArm(Lists.newArrayList(first, second, centromere));
        assertEquals(2, result.size());

        assertCombinedRegion(1, 10000, 2, CopyNumberMethod.LONG_ARM, result.get(0));
        // assertCombinedRegion(10001, 20000, 2, CopyNumberMethod.UNKNOWN, result.get(1));
    }

    static void assertCombinedRegion(int start, int end, double expectedCopyNumber, CopyNumberMethod expectedMethod,
            CombinedRegion victim)
    {
        assertEquals(expectedCopyNumber, victim.tumorCopyNumber(), EPSILON);
        assertEquals(expectedMethod, victim.copyNumberMethod());
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
    }

    private static CombinedRegion createCombinedRegion(
            final String chromosome, int start, int end, double copyNumber, double baf, SegmentSupport support)
    {
        ObservedRegion region = createDefaultFittedRegion(chromosome, start, end);
        region.setTumorCopyNumber(copyNumber);
        region.setTumorBAF(baf);
        region.setSupport(support);
        return new CombinedRegion(region);
    }

    private static CombinedRegion createCombinedRegionShortArm(int start, int end, double copyNumber, double baf, SegmentSupport support)
    {
        return createCombinedRegion(CHROMOSOME_SHORT, start, end, copyNumber, baf, support);
    }

    public static CombinedRegion createCombinedRegion(int start, int end, double copyNumber, double baf, SegmentSupport support)
    {
        return createCombinedRegion(CHROMOSOME_NORMAL, start, end, copyNumber, baf, support);
    }
}
