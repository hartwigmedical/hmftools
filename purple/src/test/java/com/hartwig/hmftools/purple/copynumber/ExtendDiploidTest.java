package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createDefaultFittedRegion;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurityAdjusterTypicalChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendDiploidTest {

    private static final int MIN_TUMOR_COUNT = 30;
    private static final int MIN_TUMOR_COUNT_AT_CENTROMERE = 50;
    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1e-10;
    private static final PurityAdjuster PURE = new PurityAdjusterTypicalChromosome(Gender.FEMALE, 1d, 1d);
    private static final ExtendDiploid PURE_VICTIM =
            new ExtendDiploid(new AlleleTolerance(PURE), MIN_TUMOR_COUNT, MIN_TUMOR_COUNT_AT_CENTROMERE);

    @Test
    public void testFavourTumorRatioCountOverLength() {
        final FittedRegion dubious1 = createDubiousRegion(1, 50000, 2, 10, 0);
        final FittedRegion dubious2 = createDubiousRegion(50001, 60000, 3, 20, 0);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(dubious1, dubious2));
        assertEquals(1, result.size());
        assertRegion(1, 60000, 3, result.get(0));
    }

    @Test
    public void testCentromereIsDubious() {
        final FittedRegion valid = createValidSomatic(1, 10000, 1, 1, SegmentSupport.TELOMERE);
        final FittedRegion dubiousCentromere = createRegion(10001, 20000, 4, 8, 0, SegmentSupport.CENTROMERE);
        final FittedRegion dubious1 = createRegion(20001, 30000, 4, 9, 0, SegmentSupport.NONE);
        final FittedRegion dubious2 = createRegion(30001, 40000, 4, 10, 0, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(valid, dubiousCentromere, dubious1, dubious2));
        assertEquals(2, result.size());
        assertRegion(1, 10000, 1, result.get(0));
        assertRegion(10001, 40000, 4, result.get(1));
    }

    @Test
    public void testDubiousRegionGetsIncludedFromLeft() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final FittedRegion dubious = createDubiousRegion(10001, 20000, 2.0, 10);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious, somaticRight));
        assertEquals(1, result.size());
        assertRegion(1, 30000, 3, result.get(0));
    }

    @Test
    public void testDubiousRegionGetsIncludedFromRight() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 30, SegmentSupport.BND);
        final FittedRegion dubious = createDubiousRegion(10001, 20000, 2.0, 10);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious, somaticRight));
        assertEquals(1, result.size());
        assertRegion(1, 30000, 3, result.get(0));
    }

    @Test
    public void testDubiousRegionGetsExcludedBecauseOfSV() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final FittedRegion dubious = createDubiousRegion(10001, 20000, 2.0, 10);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.BND);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious, somaticRight));
        assertEquals(3, result.size());
    }

    @Test
    public void testMultipleDubiousRegionsGetIncluded() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final FittedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 20);
        final FittedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 9);
        final FittedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 20);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(1, result.size());
        assertRegion(1, 30000, 3, result.get(0));
    }

    @Test
    public void testDubiousRegionGetsExcludedBecauseOfTooManyConsecutiveRatioCounts() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final FittedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 20);
        final FittedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 10);
        final FittedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 20);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(3, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 2, result.get(1));
        assertRegion(20001, 30000, 3, result.get(2));
    }

    @Test
    public void testTowardsCentromere() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final FittedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 11);
        final FittedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 11);
        final FittedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 11);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(2, result.size());
        assertRegion(1, 20000, 3, result.get(0));
        assertRegion(20001, 30000, 3, result.get(1));
    }

    @Test
    public void testTooBigTowardsCentromere() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final FittedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 21);
        final FittedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 21);
        final FittedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 21);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.CENTROMERE);

        final List<CombinedRegion> result =
                PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(3, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 2, result.get(1));
        assertRegion(20001, 30000, 3, result.get(2));
    }

    @Test
    public void testInvalidGermlineGetsIgnored() {
        final FittedRegion somatic = createFittedRegion(1, 10000, 3.0, GermlineStatus.DIPLOID, SegmentSupport.BND);
        final FittedRegion germline = createFittedRegion(10001, 20000, 4, GermlineStatus.AMPLIFICATION, SegmentSupport.NONE);
        final List<FittedRegion> regions = Lists.newArrayList(somatic, germline);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(regions);
        assertEquals(1, result.size());
        assertRegion(1, 20000, 3, result.get(0));
    }

    @Test
    public void testInvalidGermlineIsKeptWithSVSupport() {
        final FittedRegion somatic = createFittedRegion(1, 10000, 3.0, GermlineStatus.DIPLOID, SegmentSupport.BND);
        final FittedRegion germline = createFittedRegion(10001, 20000, 4.0, GermlineStatus.AMPLIFICATION, SegmentSupport.BND);
        final List<FittedRegion> regions = Lists.newArrayList(somatic, germline);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(regions);
        assertEquals(2, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 4, result.get(1));
    }

    private static void assertRegion(long start, long end, double tumorCopyNumber, @NotNull final CombinedRegion victim) {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(tumorCopyNumber, victim.tumorCopyNumber(), EPSILON);
    }

    @NotNull
    private static FittedRegion createFittedRegion(long start, long end, double tumorCopyNumber, GermlineStatus status,
            SegmentSupport support) {
        return createFittedRegion(start, end, tumorCopyNumber, 1.1, status, support);
    }

    @NotNull
    private static FittedRegion createValidSomatic(long start, long end, double copyNumber, int bafCount, SegmentSupport support) {
        return createDefaultFittedRegion(CHROMOSOME, start, end)
                .status(GermlineStatus.DIPLOID)
                .tumorCopyNumber(copyNumber)
                .refNormalisedCopyNumber(copyNumber)
                .observedNormalRatio(0.5)
                .bafCount(bafCount)
                .support(support)
                .depthWindowCount(MIN_TUMOR_COUNT_AT_CENTROMERE)
                .build();
    }

    @NotNull
    private static FittedRegion createDubiousRegion(long start, long end, double copyNumber, int ratioCount) {
        return createDubiousRegion(start, end, copyNumber, ratioCount, 1);
    }

    @NotNull
    private static FittedRegion createDubiousRegion(long start, long end, double copyNumber, int ratioCount, int bafCount) {
        return createRegion(start, end, copyNumber, ratioCount, bafCount, SegmentSupport.NONE);
    }

    @NotNull
    private static FittedRegion createRegion(long start, long end, double copyNumber, int ratioCount, int bafCount,
            SegmentSupport support) {
        return createDefaultFittedRegion(CHROMOSOME, start, end)
                .status(GermlineStatus.DIPLOID)
                .tumorCopyNumber(copyNumber)
                .refNormalisedCopyNumber(copyNumber)
                .observedNormalRatio(0.5)
                .bafCount(bafCount)
                .support(support)
                .depthWindowCount(ratioCount)
                .build();
    }

    @NotNull
    private static FittedRegion createFittedRegion(long start, long end, double tumorCopyNumber, double observedNormalRatio,
            GermlineStatus status, SegmentSupport support) {
        return createDefaultFittedRegion(CHROMOSOME, start, end)
                .status(status)
                .tumorCopyNumber(tumorCopyNumber)
                .refNormalisedCopyNumber(tumorCopyNumber)
                .observedNormalRatio(observedNormalRatio)
                .bafCount(50)
                .support(support)
                .build();
    }
}
