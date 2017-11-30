package com.hartwig.hmftools.common.purple.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ExtendDiploidTest {

    private final static int MIN_TUMOR_COUNT = 30;
    private final static String CHROMOSOME = "1";
    private final static double EPSILON = 1e-10;
    private final static PurityAdjuster PURE = new PurityAdjuster(Gender.FEMALE, 1d, 1d);
    private final static ExtendDiploid PURE_VICTIM = new ExtendDiploid(PURE, MIN_TUMOR_COUNT);

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
        final FittedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 10);
        final FittedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 9);
        final FittedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 10);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(1, result.size());
        assertRegion(1, 30000, 3, result.get(0));
    }

    @Test
    public void testDubiousRegionGetsExcludedBecauseOfTooManyConsecutiveRatioCounts() {
        final FittedRegion somaticLeft = createValidSomatic(1, 10000, 3.0, 50, SegmentSupport.BND);
        final FittedRegion dubious1 = createDubiousRegion(10001, 11000, 2.0, 10);
        final FittedRegion dubious2 = createDubiousRegion(11001, 12000, 2.0, 10);
        final FittedRegion dubious3 = createDubiousRegion(12001, 20000, 2.0, 10);
        final FittedRegion somaticRight = createValidSomatic(20001, 30000, 3.0, 40, SegmentSupport.NONE);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(Lists.newArrayList(somaticLeft, dubious1, dubious2, dubious3, somaticRight));
        assertEquals(3, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 2, result.get(1));
        assertRegion(20001, 30000, 3, result.get(2));
    }

    @Test
    public void testInvalidGermlineGetsIgnored() {
        final FittedRegion somatic = createFittedRegion(1, 10000, 3.0, ObservedRegionStatus.SOMATIC, SegmentSupport.BND);
        final FittedRegion germline = createFittedRegion(10001, 20000, 4, ObservedRegionStatus.GERMLINE_AMPLIFICATION, SegmentSupport.NONE);
        final List<FittedRegion> regions = Lists.newArrayList(somatic, germline);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(regions);
        assertEquals(1, result.size());
        assertRegion(1, 20000, 3, result.get(0));
    }

    @Test
    public void testInvalidGermlineIsKeptWithSVSupport() {
        final FittedRegion somatic = createFittedRegion(1, 10000, 3.0, ObservedRegionStatus.SOMATIC, SegmentSupport.BND);
        final FittedRegion germline =
                createFittedRegion(10001, 20000, 4.0, ObservedRegionStatus.GERMLINE_AMPLIFICATION, SegmentSupport.BND);
        final List<FittedRegion> regions = Lists.newArrayList(somatic, germline);

        final List<CombinedRegion> result = PURE_VICTIM.extendDiploid(regions);
        assertEquals(2, result.size());
        assertRegion(1, 10000, 3, result.get(0));
        assertRegion(10001, 20000, 4, result.get(1));
    }

    private void assertRegion(long start, long end, double tumorCopyNumber, @NotNull final CombinedRegion victim) {
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
        assertEquals(tumorCopyNumber, victim.tumorCopyNumber(), EPSILON);
    }

    private FittedRegion createFittedRegion(long start, long end, double tumorCopyNumber, ObservedRegionStatus status,
            SegmentSupport support) {
        return createFittedRegion(start, end, tumorCopyNumber, 1.1, status, support);
    }

    private FittedRegion createValidSomatic(long start, long end, double copyNumber, int bafCount, SegmentSupport support) {
        return PurpleDatamodelTest.createDefaultFittedRegion(CHROMOSOME, start, end)
                .status(ObservedRegionStatus.SOMATIC)
                .tumorCopyNumber(copyNumber)
                .refNormalisedCopyNumber(copyNumber)
                .observedNormalRatio(0.5)
                .bafCount(bafCount)
                .support(support)
                .observedTumorRatioCount(MIN_TUMOR_COUNT)
                .build();
    }

    private FittedRegion createDubiousRegion(long start, long end, double copyNumber, int ratioCount) {
        return PurpleDatamodelTest.createDefaultFittedRegion(CHROMOSOME, start, end)
                .status(ObservedRegionStatus.SOMATIC)
                .tumorCopyNumber(copyNumber)
                .refNormalisedCopyNumber(copyNumber)
                .observedNormalRatio(0.5)
                .bafCount(1)
                .support(SegmentSupport.NONE)
                .observedTumorRatioCount(ratioCount)
                .build();
    }

    private FittedRegion createFittedRegion(long start, long end, double tumorCopyNumber, double observedNormalRatio,
            ObservedRegionStatus status, SegmentSupport support) {
        return PurpleDatamodelTest.createDefaultFittedRegion(CHROMOSOME, start, end)
                .status(status)
                .tumorCopyNumber(tumorCopyNumber)
                .refNormalisedCopyNumber(tumorCopyNumber)
                .observedNormalRatio(observedNormalRatio)
                .bafCount(50)
                .support(support)
                .build();
    }

}
