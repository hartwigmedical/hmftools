package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class ObservedRegionStatusTest {

    @Test
    public void testStandardChromosome() {
        testStandardChromosome(Gender.MALE, "1");
        testStandardChromosome(Gender.FEMALE, "1");
    }

    @Test
    public void testSexChromosome() {
        testStandardChromosome(Gender.FEMALE, "X");
        testSpecialChromosome(Gender.FEMALE, "Y");
        testSpecialChromosome(Gender.MALE, "X");
        testSpecialChromosome(Gender.MALE, "Y");
    }

    private void testSpecialChromosome(Gender gender, String chromosome) {
        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.2));

        assertEquals(ObservedRegionStatus.SOMATIC, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.3));
        assertEquals(ObservedRegionStatus.SOMATIC, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.5));
        assertEquals(ObservedRegionStatus.SOMATIC, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.7));

        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.75));

        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.8));
        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 1));
        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 1.2));

        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 1.3));
    }


    private void testStandardChromosome(Gender gender, String chromosome) {
        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.2));

        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.3));
        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.5));
        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.7));

        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.75));

        assertEquals(ObservedRegionStatus.SOMATIC, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 0.8));
        assertEquals(ObservedRegionStatus.SOMATIC, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 1));
        assertEquals(ObservedRegionStatus.SOMATIC, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 1.2));

        assertEquals(ObservedRegionStatus.GERMLINE, ObservedRegionStatus.fromNormalRatio(gender, chromosome, 1.3));
    }

}
