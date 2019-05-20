package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;

import org.junit.Test;

public class SomaticVariantsTest {

    @Test
    public void canExtractCodingFromHGVSCodingImpactField() {
        assertEquals(927, SomaticVariants.extractCodonField("c.927+1G>A"));
        assertEquals(1799, SomaticVariants.extractCodonField("c.1799T>A"));
        assertEquals(423, SomaticVariants.extractCodonField("c.423_427delCCCTG"));
        assertEquals(8390, SomaticVariants.extractCodonField("c.8390delA"));
    }

    @Test
    public void sortCorrectlyOnCodon() {
        ReportableVariant variant1 = PatientReporterTestFactory.createTestReportableVariantBuilder().hgvsCodingImpact("c.300T>A").build();
        ReportableVariant variant2 = PatientReporterTestFactory.createTestReportableVariantBuilder().hgvsCodingImpact("c.4000T>A").build();
        ReportableVariant variant3 = PatientReporterTestFactory.createTestReportableVariantBuilder().hgvsCodingImpact("c.500T>A").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2, variant3);

        List<ReportableVariant> sortedVariants = SomaticVariants.sort(variants);

        assertEquals(variant1, sortedVariants.get(0));
        assertEquals(variant3, sortedVariants.get(1));
        assertEquals(variant2, sortedVariants.get(2));
    }

    @Test
    public void negativeMinorAllelePloidyWorks() {
        assertEquals("AAAAAAA", SomaticVariants.descriptiveBAF(7, -1));
    }

    @Test
    public void descriptiveBAFWorksAroundBoundary() {
        assertEquals("AA", SomaticVariants.descriptiveBAF(2, 2));
        assertEquals("", SomaticVariants.descriptiveBAF(-12, 0));
    }

    @Test
    public void descriptiveBAFWorksForTypicalCases() {
        assertEquals("AA", SomaticVariants.descriptiveBAF(2, 0));
        assertEquals("AB", SomaticVariants.descriptiveBAF(2, 1));
        assertEquals("AB", SomaticVariants.descriptiveBAF(2, 1.4));
        assertEquals("AA", SomaticVariants.descriptiveBAF(2, 1.6));
        assertEquals("AA", SomaticVariants.descriptiveBAF(2, 2));

        assertEquals("AAA", SomaticVariants.descriptiveBAF(3, 0));
        assertEquals("AAB", SomaticVariants.descriptiveBAF(3, 1.5));
        assertEquals("AAB", SomaticVariants.descriptiveBAF(3, 2.49));
        assertEquals("AAA", SomaticVariants.descriptiveBAF(3, 2.51));
        assertEquals("AAA", SomaticVariants.descriptiveBAF(3, 3));

        assertEquals("AAAA", SomaticVariants.descriptiveBAF(4, 0));
        assertEquals("AAAA", SomaticVariants.descriptiveBAF(4, 0.1));
        assertEquals("AABB", SomaticVariants.descriptiveBAF(4, 2));
        assertEquals("AABB", SomaticVariants.descriptiveBAF(4, 2.49));
        assertEquals("AAAB", SomaticVariants.descriptiveBAF(4, 2.51));
        assertEquals("AAAB", SomaticVariants.descriptiveBAF(4, 3.49));
        assertEquals("AAAA", SomaticVariants.descriptiveBAF(4, 3.51));
        assertEquals("AAAA", SomaticVariants.descriptiveBAF(4, 3.8));
        assertEquals("AAAA", SomaticVariants.descriptiveBAF(4, 4));

        assertEquals("A[12x]", SomaticVariants.descriptiveBAF(12, 0));
    }
}