package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.MiscTestUtils.createDefaultFittedRegion;
import static com.hartwig.hmftools.purple.MiscTestUtils.buildPurityAdjuster;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendDiploidTest
{
    private static final int MIN_TUMOR_COUNT = 30;
    private static final int MIN_TUMOR_COUNT_AT_CENTROMERE = 50;
    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1e-10;
    private static final PurityAdjuster PURE = buildPurityAdjuster(Gender.FEMALE, 1d, 1d);
    private static final ExtendDiploid PURE_VICTIM =
            new ExtendDiploid(new AlleleTolerance(PURE), MIN_TUMOR_COUNT, MIN_TUMOR_COUNT_AT_CENTROMERE);

    @Test
    public void testFavourTumorRatioCountOverLength()
    {
        final ObservedRegion dubious1 = createDubiousRegion(1, 50000, 2, 10, 0);
        final ObservedRegion dubious2 = createDubiousRegion(50001, 60000, 3, 20, 0);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(dubious1, dubious2));
        assertEquals(1, result.size());
        assertRegion(1, 60000, 3, result.get(0));
    }

    @Test
    public void testCentromereIsDubious()
    {
        final ObservedRegion valid = createValidSomatic(1, 10000, 1, 1, SegmentSupport.TELOMERE);
        final ObservedRegion dubiousCentromere = createRegion(10001, 20000, 4, 8, 0, SegmentSupport.CENTROMERE);
        final ObservedRegion dubious1 = createRegion(20001, 30000, 4, 9, 0, SegmentSupport.NONE);
        final ObservedRegion dubious2 = createRegion(30001, 40000, 4, 10, 0, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(valid, dubiousCentromere, dubious1, dubious2));
        assertEquals(2, result.size());
        assertRegion(1, 10000, 1, result.get(0));
        assertRegion(10001, 40000, 4, result.get(1));
    }

    @Test
    public void testDubiousRegionGetsIncludedFromLeft()
    {
        final ObservedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final ObservedRegion dubious = createDubiousRegion(10001, 20000, 2.0, 10);
        final ObservedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious, somaticRight));
        assertEquals(1, result.size());
        assertRegion(1, 30000, 3, result.get(0));
    }

    @Test
    public void testDubiousRegionGetsIncludedFromRight()
    {
        final ObservedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 30, SegmentSupport.BND);
        final ObservedRegion dubious = createDubiousRegion(10001, 20000, 2.0, 10);
        final ObservedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious, somaticRight));
        assertEquals(1, result.size());
        assertRegion(1, 30000, 3, result.get(0));
    }

    @Test
    public void testDubiousRegionGetsExcludedBecauseOfSV()
    {
        final ObservedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final ObservedRegion dubious = createDubiousRegion(10001, 20000, 2.0, 10);
        final ObservedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.BND);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious, somaticRight));
        assertEquals(3, result.size());
    }

    @Test
    public void testMultipleDubiousRegionsGetIncluded()
    {
        final ObservedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final ObservedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 20);
        final ObservedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 9);
        final ObservedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 20);
        final ObservedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(1, result.size());
        assertRegion(1, 30000, 3, result.get(0));
    }

    @Test
    public void testDubiousRegionGetsExcludedBecauseOfTooManyConsecutiveRatioCounts()
    {
        final ObservedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final ObservedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 20);
        final ObservedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 10);
        final ObservedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 20);
        final ObservedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(3, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 2, result.get(1));
        assertRegion(20001, 30000, 3, result.get(2));
    }

    @Test
    public void testTowardsCentromere()
    {
        final ObservedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final ObservedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 11);
        final ObservedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 11);
        final ObservedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 11);
        final ObservedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(2, result.size());
        assertRegion(1, 20000, 3, result.get(0));
        assertRegion(20001, 30000, 3, result.get(1));
    }

    @Test
    public void testTooBigTowardsCentromere()
    {
        final ObservedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final ObservedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 21);
        final ObservedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 21);
        final ObservedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 21);
        final ObservedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(3, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 2, result.get(1));
        assertRegion(20001, 30000, 3, result.get(2));
    }

    @Test
    public void testInvalidGermlineGetsIgnored()
    {
        final ObservedRegion somatic = createFittedRegion(1, 10000, 3.0, GermlineStatus.DIPLOID, SegmentSupport.BND);
        final ObservedRegion germline = createFittedRegion(10001, 20000, 4, GermlineStatus.AMPLIFICATION, SegmentSupport.NONE);
        final List<ObservedRegion> regions = Lists.newArrayList(somatic, germline);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(regions);
        assertEquals(1, result.size());
        assertRegion(1, 20000, 3, result.get(0));
    }

    @Test
    public void testInvalidGermlineIsKeptWithSVSupport()
    {
        final ObservedRegion somatic = createFittedRegion(1, 10000, 3.0, GermlineStatus.DIPLOID, SegmentSupport.BND);
        final ObservedRegion germline = createFittedRegion(10001, 20000, 4.0, GermlineStatus.AMPLIFICATION, SegmentSupport.BND);
        final List<ObservedRegion> regions = Lists.newArrayList(somatic, germline);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(regions);
        assertEquals(2, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 4, result.get(1));
    }

    private static void assertRegion(int start, int end, double tumorCopyNumber, @NotNull final CombinedRegion victim)
    {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(tumorCopyNumber, victim.tumorCopyNumber(), EPSILON);
    }

    @NotNull
    private static ObservedRegion createFittedRegion(int start, int end, double tumorCopyNumber, GermlineStatus status,
            SegmentSupport support)
    {
        return createFittedRegion(start, end, tumorCopyNumber, 1.1, status, support);
    }

    @NotNull
    private static ObservedRegion createValidSomatic(int start, int end, double copyNumber, int bafCount, SegmentSupport support)
    {
        ObservedRegion region = createDefaultFittedRegion(CHROMOSOME, start, end);
        region.setGermlineStatus(GermlineStatus.DIPLOID);
        region.setTumorCopyNumber(copyNumber);
        region.setRefNormalisedCopyNumber(copyNumber);
        region.setObservedNormalRatio(0.5);
        region.setBafCount(bafCount);
        region.setSupport(support);
        region.setDepthWindowCount(MIN_TUMOR_COUNT_AT_CENTROMERE);
        return region;
    }

    @NotNull
    private static ObservedRegion createDubiousRegion(int start, int end, double copyNumber, int ratioCount)
    {
        return createDubiousRegion(start, end, copyNumber, ratioCount, 1);
    }

    @NotNull
    private static ObservedRegion createDubiousRegion(int start, int end, double copyNumber, int ratioCount, int bafCount)
    {
        return createRegion(start, end, copyNumber, ratioCount, bafCount, SegmentSupport.NONE);
    }

    @NotNull
    private static ObservedRegion createRegion(int start, int end, double copyNumber, int ratioCount, int bafCount,
            SegmentSupport support)
    {
        ObservedRegion region = createDefaultFittedRegion(CHROMOSOME, start, end);
        region.setGermlineStatus(GermlineStatus.DIPLOID);
        region.setTumorCopyNumber(copyNumber);
        region.setRefNormalisedCopyNumber(copyNumber);
        region.setObservedNormalRatio(0.5);
        region.setBafCount(bafCount);
        region.setSupport(support);
        region.setDepthWindowCount(ratioCount);
        return region;
    }

    @NotNull
    private static ObservedRegion createFittedRegion(int start, int end, double tumorCopyNumber, double observedNormalRatio,
            GermlineStatus status, SegmentSupport support)
    {
        ObservedRegion region = createDefaultFittedRegion(CHROMOSOME, start, end);
        region.setGermlineStatus(status);
        region.setTumorCopyNumber(tumorCopyNumber);
        region.setRefNormalisedCopyNumber(tumorCopyNumber);
        region.setObservedNormalRatio(observedNormalRatio);
        region.setBafCount(50);
        region.setSupport(support);
        return region;
    }
}
