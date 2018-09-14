package com.hartwig.hmftools.patientreporter.report.data;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class VariantDataSourceTest {

    @Test
    public void negativeMinorAllelePloidyWorks() {
        assertEquals("AAAAAAA", VariantDataSource.descriptiveBAF(7, -1));
    }

    @Test
    public void descriptiveBAFWorksAroundBoundary() {
        assertEquals("AA", VariantDataSource.descriptiveBAF(2, 2));
        assertEquals("", VariantDataSource.descriptiveBAF(-12, 0));
    }

    @Test
    public void descriptiveBAFWorksForTypicalCases() {
        assertEquals("AA", VariantDataSource.descriptiveBAF(2, 0));
        assertEquals("AB", VariantDataSource.descriptiveBAF(2, 1));
        assertEquals("AB", VariantDataSource.descriptiveBAF(2, 1.4));
        assertEquals("AA", VariantDataSource.descriptiveBAF(2, 1.6));
        assertEquals("AA", VariantDataSource.descriptiveBAF(2, 2));

        assertEquals("AAA", VariantDataSource.descriptiveBAF(3, 0));
        assertEquals("AAB", VariantDataSource.descriptiveBAF(3, 1.5));
        assertEquals("AAB", VariantDataSource.descriptiveBAF(3, 2.49));
        assertEquals("AAA", VariantDataSource.descriptiveBAF(3, 2.51));
        assertEquals("AAA", VariantDataSource.descriptiveBAF(3, 3));

        assertEquals("AAAA", VariantDataSource.descriptiveBAF(4, 0));
        assertEquals("AAAA", VariantDataSource.descriptiveBAF(4, 0.1));
        assertEquals("AABB", VariantDataSource.descriptiveBAF(4, 2));
        assertEquals("AABB", VariantDataSource.descriptiveBAF(4, 2.49));
        assertEquals("AAAB", VariantDataSource.descriptiveBAF(4, 2.51));
        assertEquals("AAAB", VariantDataSource.descriptiveBAF(4, 3.49));
        assertEquals("AAAA", VariantDataSource.descriptiveBAF(4, 3.51));
        assertEquals("AAAA", VariantDataSource.descriptiveBAF(4, 3.8));
        assertEquals("AAAA", VariantDataSource.descriptiveBAF(4, 4));

        assertEquals("A[12x]", VariantDataSource.descriptiveBAF(12, 0));
    }
}