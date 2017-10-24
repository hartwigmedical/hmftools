package com.hartwig.hmftools.common.purple.segment;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.junit.Test;

public class SegmentStatusTest {

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
        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.2));

        assertEquals(SegmentStatus.SOMATIC, SegmentStatus.fromNormalRatio(gender, chromosome, 0.3));
        assertEquals(SegmentStatus.SOMATIC, SegmentStatus.fromNormalRatio(gender, chromosome, 0.5));
        assertEquals(SegmentStatus.SOMATIC, SegmentStatus.fromNormalRatio(gender, chromosome, 0.7));

        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.75));

        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.8));
        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 1));
        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 1.2));

        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 1.3));
    }


    private void testStandardChromosome(Gender gender, String chromosome) {
        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.2));

        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.3));
        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.5));
        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.7));

        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 0.75));

        assertEquals(SegmentStatus.SOMATIC, SegmentStatus.fromNormalRatio(gender, chromosome, 0.8));
        assertEquals(SegmentStatus.SOMATIC, SegmentStatus.fromNormalRatio(gender, chromosome, 1));
        assertEquals(SegmentStatus.SOMATIC, SegmentStatus.fromNormalRatio(gender, chromosome, 1.2));

        assertEquals(SegmentStatus.GERMLINE, SegmentStatus.fromNormalRatio(gender, chromosome, 1.3));
    }

}
