package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.BAFUtils.minAlleleCount;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BAFUtilsTest {

    private static final double EPSILON = 1e-10;

    private static final double NORMAL_BAF = 0.535;
    private static final BAFUtils BAF_UTILS = new BAFUtils(NORMAL_BAF);

    @Test
    public void testMinBetaAllele() {
        assertEquals(0, minAlleleCount(-28));

        assertEquals(1, minAlleleCount(1));
        assertEquals(1, minAlleleCount(2));
        assertEquals(2, minAlleleCount(3));
        assertEquals(2, minAlleleCount(4));
        assertEquals(3, minAlleleCount(5));
        assertEquals(3, minAlleleCount(6));
        assertEquals(4, minAlleleCount(7));
    }

    @Test
    public void testDiploidModelBaf() {
        assertEquals(NORMAL_BAF, BAF_UTILS.modelBAF(1, 2, 1), EPSILON);
        assertEquals(NORMAL_BAF, BAF_UTILS.modelBAF(0.5, 4, 2), EPSILON);
        assertEquals(NORMAL_BAF, BAF_UTILS.modelBAF(0, 6, 3), EPSILON);
    }

    @Test
    public void testZeroPloidyModelBaf() {
        assertEquals(NORMAL_BAF, BAF_UTILS.modelBAF(0, 0, 0), EPSILON);
        assertEquals(NORMAL_BAF, BAF_UTILS.modelBAF(0.5, 0, 0), EPSILON);
        assertEquals(NORMAL_BAF, BAF_UTILS.modelBAF(1, 0, 0), EPSILON);
    }

    @Test
    public void testIsClonal() {
        assertTrue(BAFUtils.isClonal(0.75));
        assertTrue(BAFUtils.isClonal(0.9));
        assertTrue(BAFUtils.isClonal(1.0));
        assertTrue(BAFUtils.isClonal(1.1));
        assertTrue(BAFUtils.isClonal(1.25));

        assertFalse(BAFUtils.isClonal(1.26));
    }

    @Test
    public void testBaf() {
        assertBaf(1, 0.8, 0.5, 0.6);
        assertBaf(1, 0.8, 1, 0.6);

        assertBaf(0.625, 0.8, 2, 0.6);
        assertBafMaleSex(1, 0.8, 2, 0.6);
    }

    @Test
    public void testBafAdjustment() {
        assertBaf(0.60, 0.08, 5.1, 0.5);
        assertBaf(0.60, 0.08, 5.2, 0.5);
        assertBaf(0.60, 0.08, 5.25, 0.5);
        assertBaf(0.50, 0.08, 5.5, 0.5);
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

    @Test
    public void testPurityAdjustedVAF() {
        testPurityAdjustedVAF(Gender.MALE, "X", 0.3, 2, 2);
        testPurityAdjustedVAF(Gender.MALE, "X", 0.3, 2, 1);
        testPurityAdjustedVAF(Gender.FEMALE, "X", 0.3, 2, 1);
        testPurityAdjustedVAF(Gender.FEMALE, "X", 0.4, 3, 1);
        testPurityAdjustedVAF(Gender.FEMALE, "X", 0.4, 4, 1);
        testPurityAdjustedVAF(Gender.FEMALE, "X", 0.4, 4, 2);
    }

    @Test
    public void testObservedVAF() {
        final double male = observedVAF(Gender.MALE, "1", 0.5, 2, 2);
        final double maleX = observedVAF(Gender.MALE, "X", 0.5, 2, 1);
        final double femaleX = observedVAF(Gender.FEMALE, "X", 0.5, 2, 1);

        assertEquals(1/3d, maleX, EPSILON);
        assertEquals(1/4d, femaleX, EPSILON);

        assertEquals(1/2d, male, EPSILON);
    }

    private static double expectedFrequency(final int ploidy, final int alleleCount) {
        return 1d * alleleCount / ploidy;
    }

    private static double observedVAF(@NotNull final Gender gender, @NotNull final String chromosome, final double purity, final int ploidy,
            final int alleleCount) {
        double expectedPurityAdjustedVAF = 1d * alleleCount / ploidy;
        int normalCount = HumanChromosome.fromString(chromosome).isDiploid(gender) ? 2 : 1;
        return expectedPurityAdjustedVAF * ploidy * purity / (normalCount * (1 - purity) + ploidy * purity);
    }

    private void testPurityAdjustedVAF(@NotNull final Gender gender, @NotNull final String chromosome, final double purity, final int ploidy,
            final int alleleCount) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, 1);
        double expected = expectedFrequency(ploidy, alleleCount);
        double observed = observedVAF(gender, chromosome, purity, ploidy, alleleCount);
        assertEquals(expected, purityAdjuster.purityAdjustedVAF(chromosome, ploidy, observed), EPSILON);
    }

    private void testPurityAdjustedBaf(final double purity, final int ploidy, final int alleleCount) {
        PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.FEMALE, purity, 1);
        double expectedPurityAdjustedBAF = 1d * alleleCount / ploidy;
        double observedBAF = modelBAF(purity, ploidy, alleleCount);
        assertEquals(expectedPurityAdjustedBAF, BAF_UTILS.purityAdjustedBAF(purityAdjuster,"1", ploidy, observedBAF), EPSILON);
    }

    private void assertBaf(final double expectedBaf, final Gender gender, final String chromsome, final double purity,
            final double copyNumber, final double observedFrequency) {
        PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, 1);
        assertEquals(expectedBaf, BAF_UTILS.purityAdjustedBAF(purityAdjuster, chromsome, copyNumber, observedFrequency), EPSILON);
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


}
