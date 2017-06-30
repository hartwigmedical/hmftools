package com.hartwig.hmftools.common.purity;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class PurityAdjusterTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testBaf() {
        assertBaf( 0.375, 0.8, 2, 0.4 );
        assertBafMaleSex( 0.5, 0.8, 2, 0.4 );
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

    private void assertBaf(final double expectedBaf, final Gender gender, final String chromsome, final double purity, final double copyNumber, final double observedFrequency) {
        PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, 1);
        assertEquals(expectedBaf, purityAdjuster.purityAdjustedBAF(chromsome, copyNumber, observedFrequency), EPSILON);
    }

    @Test
    public void testPurityAdjustedCopyNumber() {
        assertPurityAdjustment(0, 0.85, 1.04, 0);
        assertPurityAdjustment(1, 0.85, 1.0, 0.575);
    }

    private void assertPurityAdjustment(final double expectedAdjustedCopyNumber, final double purity, final double normFactor, final double ratio) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.MALE, purity, normFactor);
        assertEquals(expectedAdjustedCopyNumber, purityAdjuster.purityAdjustedCopyNumber("1", ratio), EPSILON);
    }

    private PurityAdjuster adjuster(final double purity, final double normFactor) {
        return new PurityAdjuster(Gender.MALE, purity, normFactor);
    }

}
