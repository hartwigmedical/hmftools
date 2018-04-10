package com.hartwig.hmftools.common.purple.gender;

import static com.hartwig.hmftools.common.purple.gender.Gender.FEMALE;
import static com.hartwig.hmftools.common.purple.gender.Gender.MALE;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.cobalt.ImmutableReferenceRatioStatistics;
import com.hartwig.hmftools.common.cobalt.ReferenceRatioStatistics;

import org.junit.Test;

public class GenderTest {

    @Test
    public void testGender() {

        assertGender(Gender.MALE_KLINEFELTER, 0.98, 0.48, 9842);

        assertGender(FEMALE, 0.6872, 0.0049, 1304);
        assertGender(FEMALE, 0.74, 0.01, 947);

        assertGender(MALE, 0.50, 0.08, 9837);
        assertGender(MALE, 0.51, 0.47, 9842);
        assertGender(MALE, 0.51, 1.02, 9842);
    }

    private void assertGender(Gender gender, double xMedian, double yMedian, int yCount) {
        final ReferenceRatioStatistics statistic = create(xMedian, yMedian, yCount);
        assertEquals(gender, Gender.fromCobalt(statistic));
    }

    private ReferenceRatioStatistics create(double xMedian, double yMedian, int yCount) {
        return ImmutableReferenceRatioStatistics.builder().xCount(1).xMedian(xMedian).yMedian(yMedian).yCount(yCount).build();
    }

}