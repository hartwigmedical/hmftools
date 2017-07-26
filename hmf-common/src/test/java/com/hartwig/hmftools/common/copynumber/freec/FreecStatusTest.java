package com.hartwig.hmftools.common.copynumber.freec;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class FreecStatusTest {

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
        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.2));

        assertEquals(FreecStatus.SOMATIC, FreecStatus.fromNormalRatio(gender, chromosome, 0.3));
        assertEquals(FreecStatus.SOMATIC, FreecStatus.fromNormalRatio(gender, chromosome, 0.5));
        assertEquals(FreecStatus.SOMATIC, FreecStatus.fromNormalRatio(gender, chromosome, 0.7));

        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.75));

        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.8));
        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 1));
        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 1.2));

        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 1.3));
    }


    private void testStandardChromosome(Gender gender, String chromosome) {
        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.2));

        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.3));
        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.5));
        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.7));

        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 0.75));

        assertEquals(FreecStatus.SOMATIC, FreecStatus.fromNormalRatio(gender, chromosome, 0.8));
        assertEquals(FreecStatus.SOMATIC, FreecStatus.fromNormalRatio(gender, chromosome, 1));
        assertEquals(FreecStatus.SOMATIC, FreecStatus.fromNormalRatio(gender, chromosome, 1.2));

        assertEquals(FreecStatus.GERMLINE, FreecStatus.fromNormalRatio(gender, chromosome, 1.3));
    }

}
