package com.hartwig.hmftools.common.purple.copynumber;

import static com.hartwig.hmftools.common.numeric.Doubles.greaterThan;
import static com.hartwig.hmftools.common.numeric.Doubles.lessThan;
import static com.hartwig.hmftools.common.purple.copynumber.HighConfidenceSmoothedRegions.MIN_RATIO_ONLY_TUMOR_RATIO_COUNT;
import static com.hartwig.hmftools.common.purple.copynumber.HighConfidenceSmoothedRegions.allowedBAFDeviation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.junit.Test;

public class HighConfidenceSmoothedRegionsTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void bafAllowances() {
        assertTrue(lessThan(allowedBAFDeviation(1), 0.44));
        assertTrue(greaterThan(allowedBAFDeviation(7), 0.072));
        assertTrue(greaterThan(allowedBAFDeviation(14), 0.052));
        assertTrue(greaterThan(allowedBAFDeviation(46), 0.03));
    }

    @Test
    public void beforeFirstHighConfidenceRegion() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2));
        final List<FittedRegion> copyNumbers = Lists.newArrayList(createFittedCopyNumber(1, 300, 1),
                createFittedCopyNumber(301, 1000, 2),
                createFittedCopyNumber(1001, 2000, 2));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(2, results.size());

        assertRegion(results.get(0), 1, 300, 1);
        assertRegion(results.get(1), 301, 2000, 2);
    }

    @Test
    public void betweenHighConfidenceRegions() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2), createRegion(3001, 4000, 3));
        final List<FittedRegion> copyNumbers = Lists.newArrayList(createFittedCopyNumber(1001, 2000, 2),
                createFittedCopyNumber(2001, 2200, 2),
                createFittedCopyNumber(2201, 2500, 1),
                createFittedCopyNumber(2501, 3000, 3),
                createFittedCopyNumber(3001, 4000, 3));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(3, results.size());

        assertRegion(results.get(0), 1001, 2200, 2);
        assertRegion(results.get(1), 2201, 2500, 1);
        assertRegion(results.get(2), 2501, 4000, 3);
    }

    @Test
    public void inHighConfidenceRegions() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2));

        final List<FittedRegion> copyNumbers = Lists.newArrayList(createFittedCopyNumber(1001, 1110, 2),
                createFittedCopyNumber(1111, 1220, 1),
                createFittedCopyNumber(1221, 1330, 1),
                createFittedCopyNumber(1331, 1440, 3),
                createFittedCopyNumber(1441, 2000, 2));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(4, results.size());

        assertRegion(results.get(0), 1001, 1110, 2);
        assertRegion(results.get(1), 1111, 1330, 1);
        assertRegion(results.get(2), 1331, 1440, 3);
        assertRegion(results.get(3), 1441, 2000, 2);
    }

    @Test
    public void finalPass() {
        final double baseCopyNumber = 2;

        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(createRegion(1001, 2000, baseCopyNumber));

        final List<FittedRegion> copyNumbers = Lists.newArrayList(createFittedCopyNumber(1001, 1110, baseCopyNumber),
                createFittedCopyNumber(1111, 1220, baseCopyNumber + CopyNumberDeviation.MIN_COPY_NUMBER_TOLERANCE + 0.1),
                createFittedCopyNumber(1221, 2000, baseCopyNumber + CopyNumberDeviation.MIN_COPY_NUMBER_TOLERANCE - 0.1));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(1, results.size());

        assertRegion(results.get(0), 1001, 2000, 2.2);
    }

    @Test
    public void afterLastHighConfidenceRegion() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2));
        final List<FittedRegion> copyNumbers = Lists.newArrayList(createFittedCopyNumber(1001, 2000, 2),
                createFittedCopyNumber(2001, 3000, 2),
                createFittedCopyNumber(3001, 5000, 5));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(2, results.size());

        assertRegion(results.get(0), 1001, 3000, 2);
        assertRegion(results.get(1), 3001, 5000, 5);
    }

    @Test
    public void testFavourStructuralSupportInGermlineEarly() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2), createRegion(3001, 4000, 3));

        final List<FittedRegion> copyNumbers = Lists.newArrayList(createFittedCopyNumber(1001, 2000, 2),
                createGermlineRegion(2001, 2500, 2.5, StructuralVariantSupport.DEL),
                createGermlineRegion(2501, 3000, 2.5, StructuralVariantSupport.NONE),
                createFittedCopyNumber(3001, 4000, 3));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(2, results.size());

        assertRegion(results.get(0), 1001, 2000, 2);
        assertRegion(results.get(1), 2001, 4000, 3);
    }

    @Test
    public void testFavourStructuralSupportInGermlineLate() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(createRegion(1001, 2000, 2), createRegion(3001, 4000, 3));

        final List<FittedRegion> copyNumbers = Lists.newArrayList(createFittedCopyNumber(1001, 2000, 2),
                createGermlineRegion(2001, 2500, 2.5, StructuralVariantSupport.DEL),
                createFittedCopyNumber(2501, 3000, 3),
                createFittedCopyNumber(3001, 4000, 3));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(2, results.size());

        assertRegion(results.get(0), 1001, 2000, 2);
        assertRegion(results.get(1), 2001, 4000, 3);
    }

    @Test
    public void mergeLowTumorCountSegmentStart() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(
                createRegion(1001, 2000, 1),
                createRegion(3001, 4000, 3));

        final List<FittedRegion> copyNumbers = Lists.newArrayList(
                    createFittedCopyNumber(1001, 2000, 1, MIN_RATIO_ONLY_TUMOR_RATIO_COUNT - 1),
                    createFittedCopyNumber(2001, 3000, 2),
                    createFittedCopyNumber(3001, 4000, 3),
                    createFittedCopyNumber(4001, 5000, 4));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(3, results.size());

        assertRegion(results.get(0), 1001, 3000, 2);
        assertRegion(results.get(1), 3001, 4000, 3);
        assertRegion(results.get(2), 4001, 5000, 4);
    }

    @Test
    public void mergeLowTumorCountSegment() {
        final List<PurpleCopyNumber> broadRegions = Lists.newArrayList(
                createRegion(1001, 2000, 1),
                createRegion(3001, 4000, 3));

        final List<FittedRegion> copyNumbers = Lists.newArrayList(
                createFittedCopyNumber(1001, 2000, 1),
                createFittedCopyNumber(2001, 3000, 2, MIN_RATIO_ONLY_TUMOR_RATIO_COUNT - 1),
                createFittedCopyNumber(3001, 4000, 3),
                createFittedCopyNumber(4001, 5000, 4));

        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, 1, 1);
        final List<FittedRegion> results = new HighConfidenceSmoothedRegions(purityAdjuster, broadRegions, copyNumbers).smoothedRegions();
        assertEquals(3, results.size());

        assertRegion(results.get(0), 1001, 3000, 1);
        assertRegion(results.get(1), 3001, 4000, 3);
        assertRegion(results.get(2), 4001, 5000, 4);
    }

    private void assertRegion(FittedRegion victim, long expectedStart, long expectedEnd, double expectedCopyNumber) {
        assertEquals(expectedStart, victim.start());
        assertEquals(expectedEnd, victim.end());
        assertEquals(expectedCopyNumber, victim.tumorCopyNumber(), EPSILON);
    }

    private FittedRegion createFittedCopyNumber(long start, long end, double copyNumber, int tumorCount) {
        return PurpleDatamodelTest.createDefaultFittedRegion("1", start, end)
                .bafCount(CopyNumberDeviation.MAX_BAF_COUNT)
                .observedBAF(0.5)
                .tumorCopyNumber(copyNumber)
                .refNormalisedCopyNumber(copyNumber)
                .observedTumorRatioCount(tumorCount)
                .status(ObservedRegionStatus.SOMATIC)
                .build();
    }

    private FittedRegion createFittedCopyNumber(long start, long end, double copyNumber) {
        return PurpleDatamodelTest.createDefaultFittedRegion("1", start, end)
                .bafCount(CopyNumberDeviation.MAX_BAF_COUNT)
                .observedBAF(0.5)
                .tumorCopyNumber(copyNumber)
                .refNormalisedCopyNumber(copyNumber)
                .observedTumorRatioCount(MIN_RATIO_ONLY_TUMOR_RATIO_COUNT)
                .build();
    }

    private FittedRegion createGermlineRegion(long start, long end, double copyNumber, StructuralVariantSupport sv) {
        return PurpleDatamodelTest.createDefaultFittedRegion("1", start, end)
                .bafCount(CopyNumberDeviation.MAX_BAF_COUNT)
                .status(ObservedRegionStatus.GERMLINE)
                .structuralVariantSupport(sv)
                .observedBAF(0.5)
                .tumorCopyNumber(copyNumber)
                .refNormalisedCopyNumber(copyNumber)
                .observedTumorRatioCount(MIN_RATIO_ONLY_TUMOR_RATIO_COUNT)
                .build();
    }

    private PurpleCopyNumber createRegion(long start, long end, double copyNumber) {
        return PurpleDatamodelTest.createCopyNumber("1", start, end, copyNumber)
                .bafCount(CopyNumberDeviation.MAX_BAF_COUNT)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5)
                .build();
    }
}
