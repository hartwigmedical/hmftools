package com.hartwig.hmftools.sage.coverage;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GeneCoverageTest {
    @Test
    public void testBucket() {
        for (int depth = 0; depth < 30; depth++) {
            assertEquals(depth, GeneCoverage.bucket(depth));
        }

        assertEquals(30, GeneCoverage.bucket(30));
        assertEquals(30, GeneCoverage.bucket(31));
        assertEquals(30, GeneCoverage.bucket(32));
        assertEquals(30, GeneCoverage.bucket(33));
        assertEquals(30, GeneCoverage.bucket(34));
        assertEquals(30, GeneCoverage.bucket(35));
        assertEquals(30, GeneCoverage.bucket(36));
        assertEquals(30, GeneCoverage.bucket(37));
        assertEquals(30, GeneCoverage.bucket(38));
        assertEquals(30, GeneCoverage.bucket(39));

        assertEquals(31, GeneCoverage.bucket(40));
        assertEquals(31, GeneCoverage.bucket(49));
        assertEquals(32, GeneCoverage.bucket(50));
        assertEquals(33, GeneCoverage.bucket(60));
        assertEquals(34, GeneCoverage.bucket(70));
        assertEquals(35, GeneCoverage.bucket(80));
        assertEquals(36, GeneCoverage.bucket(90));
        assertEquals(37, GeneCoverage.bucket(100));
        assertEquals(37, GeneCoverage.bucket(101));
        assertEquals(37, GeneCoverage.bucket(150));
    }

}
