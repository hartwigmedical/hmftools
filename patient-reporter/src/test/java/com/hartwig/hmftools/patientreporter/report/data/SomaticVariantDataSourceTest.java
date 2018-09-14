package com.hartwig.hmftools.patientreporter.report.data;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SomaticVariantDataSourceTest {

    @Test
    public void negativeMinorAllelePloidyWorks() {
        assertEquals("AAAAAAA", SomaticVariantDataSource.descriptiveBAF(7, -1));
    }

    @Test
    public void descriptiveBAFWorksAroundBoundary() {
        assertEquals("AA", SomaticVariantDataSource.descriptiveBAF(2, 2));
        assertEquals("", SomaticVariantDataSource.descriptiveBAF(-12, 0));
    }

    @Test
    public void descriptiveBAFWorksForTypicalCases() {
        assertEquals("AA", SomaticVariantDataSource.descriptiveBAF(2, 0));
        assertEquals("AB", SomaticVariantDataSource.descriptiveBAF(2, 1));
        assertEquals("AB", SomaticVariantDataSource.descriptiveBAF(2, 1.4));
        assertEquals("AA", SomaticVariantDataSource.descriptiveBAF(2, 1.6));
        assertEquals("AA", SomaticVariantDataSource.descriptiveBAF(2, 2));

        assertEquals("AAA", SomaticVariantDataSource.descriptiveBAF(3, 0));
        assertEquals("AAB", SomaticVariantDataSource.descriptiveBAF(3, 1.5));
        assertEquals("AAB", SomaticVariantDataSource.descriptiveBAF(3, 2.49));
        assertEquals("AAA", SomaticVariantDataSource.descriptiveBAF(3, 2.51));
        assertEquals("AAA", SomaticVariantDataSource.descriptiveBAF(3, 3));

        assertEquals("AAAA", SomaticVariantDataSource.descriptiveBAF(4, 0));
        assertEquals("AAAA", SomaticVariantDataSource.descriptiveBAF(4, 0.1));
        assertEquals("AABB", SomaticVariantDataSource.descriptiveBAF(4, 2));
        assertEquals("AABB", SomaticVariantDataSource.descriptiveBAF(4, 2.49));
        assertEquals("AAAB", SomaticVariantDataSource.descriptiveBAF(4, 2.51));
        assertEquals("AAAB", SomaticVariantDataSource.descriptiveBAF(4, 3.49));
        assertEquals("AAAA", SomaticVariantDataSource.descriptiveBAF(4, 3.51));
        assertEquals("AAAA", SomaticVariantDataSource.descriptiveBAF(4, 3.8));
        assertEquals("AAAA", SomaticVariantDataSource.descriptiveBAF(4, 4));

        assertEquals("A[12x]", SomaticVariantDataSource.descriptiveBAF(12, 0));
    }
}