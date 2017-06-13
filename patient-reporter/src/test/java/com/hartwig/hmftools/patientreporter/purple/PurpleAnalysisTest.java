package com.hartwig.hmftools.patientreporter.purple;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.junit.Test;

public class PurpleAnalysisTest {

    @Test
    public void testDiploid() {
        testAdjustedCopyNumber(2, 0.1, 0);
        testAdjustedCopyNumber(2, 0.9, 1);
        testAdjustedCopyNumber(2, 2.1, 2);
        testAdjustedCopyNumber(2, 2.9, 3);
    }

    @Test
    public void testAnuploidy() {
        testAdjustedCopyNumber(0.7, 0.1, 0);
        testAdjustedCopyNumber(0.7, 0.9, 1);
        testAdjustedCopyNumber(0.7, 2.1, 6);
        testAdjustedCopyNumber(0.7, 2.9, 8);

        testAdjustedCopyNumber(1.5, 0.1, 0);
        testAdjustedCopyNumber(1.5, 0.9, 1);
        testAdjustedCopyNumber(1.5, 2.1, 3);
        testAdjustedCopyNumber(1.5, 2.9, 4);
    }

    private void testAdjustedCopyNumber(double ploidy, double copyNumber, int expectedCopyNumber) {
        final CopyNumber adjustedCopyNumber = PurpleAnalysis.ploidyAdjusted(ploidy, create(copyNumber));
        assertCopyNumber(expectedCopyNumber, adjustedCopyNumber);
    }

    private void assertCopyNumber(int expected, CopyNumber copyNumber) {
        assertEquals(expected, copyNumber.value());
    }



    private static PurpleCopyNumber create(double averageTumorCopyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome("1")
                .start(1)
                .end(100)
                .averageTumorCopyNumber(averageTumorCopyNumber)
                .averageActualBAF(0.5)
                .averageObservedBAF(0.5)
                .bafCount(20)
                .build();
    }

}
