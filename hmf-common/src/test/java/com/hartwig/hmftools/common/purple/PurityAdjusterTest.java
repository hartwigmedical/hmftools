package com.hartwig.hmftools.common.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class PurityAdjusterTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testBaf() {
        assertBaf(0, 0.8, 0.5, 0.6);
        assertBaf(1, 0.8, 1, 0.6);

        assertBaf(0.625, 0.8, 2, 0.6);
        assertBafMaleSex(0, 0.8, 2, 0.6);
    }

    @Test
    public void testBafAdjustment() {
        assertBaf(0.60, 0.08, 5.1, 0.5);
        assertBaf(0.60, 0.08, 5.2, 0.5);
        assertBaf(0.60, 0.08, 5.25, 0.5);
        assertBaf(0.50, 0.08, 5.5, 0.5);
    }

    @Test
    public void testIsClonal() {
        assertTrue(PurityAdjuster.isClonal(0.75));
        assertTrue(PurityAdjuster.isClonal(0.9));
        assertTrue(PurityAdjuster.isClonal(1.0));
        assertTrue(PurityAdjuster.isClonal(1.1));
        assertTrue(PurityAdjuster.isClonal(1.25));

        assertFalse(PurityAdjuster.isClonal(1.26));
    }

    @Test
    public void testPurityAdjustedBaf() {
        testPurityAdjustedBaf(0.1);
        testPurityAdjustedBaf(0.2);
        testPurityAdjustedBaf(0.3);
        testPurityAdjustedBaf(0.4);
        testPurityAdjustedBaf(0.5);
        testPurityAdjustedBaf(0.6);
        testPurityAdjustedBaf(0.7);
        testPurityAdjustedBaf(0.8);
        testPurityAdjustedBaf(0.9);
    }

    private void testPurityAdjustedBaf(double purity) {
        testPurityAdjustedBaf(purity, 2, 1);
        testPurityAdjustedBaf(purity, 2, 2);
        testPurityAdjustedBaf(purity, 3, 2);
        testPurityAdjustedBaf(purity, 3, 3);
        testPurityAdjustedBaf(purity, 4, 2);
        testPurityAdjustedBaf(purity, 4, 3);
        testPurityAdjustedBaf(purity, 4, 4);
        testPurityAdjustedBaf(purity, 5, 3);
        testPurityAdjustedBaf(purity, 5, 4);
    }

    private void testPurityAdjustedBaf(final double purity, final int ploidy, final int alleleCount) {
        PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.FEMALE, purity, 1);
        double expectedPurityAdjustedBAF = 1d * alleleCount / ploidy;
        double observedBAF = modelBAF(purity, ploidy, alleleCount);
        assertEquals(expectedPurityAdjustedBAF, purityAdjuster.purityAdjustedBAF("1", ploidy, observedBAF), EPSILON);
    }

    private static double modelBAF(final double purity, final int ploidy, final int alleleCount) {
        assert (alleleCount >= ploidy / 2d);
        return (1 + purity * (alleleCount - 1)) / (2 + purity * (ploidy - 2));
    }

    private void assertBaf(final double expectedBaf, final double purity, final double copyNumber, final double observedFrequency) {
        assertBaf(expectedBaf, Gender.FEMALE, "X", purity, copyNumber, observedFrequency);
        assertBaf(expectedBaf, Gender.FEMALE, "1", purity, copyNumber, observedFrequency);
        assertBaf(expectedBaf, Gender.MALE, "1", purity, copyNumber, observedFrequency);
    }

    private void assertBafMaleSex(final double expectedBaf, final double purity, final double copyNumber, final double observedFrequency) {
        assertBaf(expectedBaf, Gender.MALE, "X", purity, copyNumber, observedFrequency);
        assertBaf(expectedBaf, Gender.MALE, "Y", purity, copyNumber, observedFrequency);
    }

    private void assertBaf(final double expectedBaf, final Gender gender, final String chromsome, final double purity,
            final double copyNumber, final double observedFrequency) {
        PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, 1);
        assertEquals(expectedBaf, purityAdjuster.purityAdjustedBAF(chromsome, copyNumber, observedFrequency), EPSILON);
    }

    @Test
    public void testPurityAdjustedCopyNumber() {
        assertPurityAdjustment(0, 0.85, 1.04, 0);
        assertPurityAdjustment(1, 0.85, 1.0, 0.575);
    }

    private void assertPurityAdjustment(final double expectedAdjustedCopyNumber, final double purity, final double normFactor,
            final double ratio) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, purity, normFactor);
        assertEquals(expectedAdjustedCopyNumber, purityAdjuster.purityAdjustedCopyNumber("1", ratio), EPSILON);
    }

    private PurityAdjuster adjuster(final double purity, final double normFactor) {
        return new PurityAdjuster(Gender.MALE, purity, normFactor);
    }

}
