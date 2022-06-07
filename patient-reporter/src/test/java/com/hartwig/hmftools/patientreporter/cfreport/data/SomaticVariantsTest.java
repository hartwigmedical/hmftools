package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;

import org.junit.Test;

public class SomaticVariantsTest {

    @Test
    public void canExtractCodingFromHGVSCodingImpactField() {
        assertEquals(927, SomaticVariants.extractCodonField("c.927+1G>A"));
        assertEquals(1799, SomaticVariants.extractCodonField("c.1799T>A"));
        assertEquals(423, SomaticVariants.extractCodonField("c.423_427delCCCTG"));
        assertEquals(8390, SomaticVariants.extractCodonField("c.8390delA"));
        assertEquals(-124, SomaticVariants.extractCodonField("c.-124C>T"));
    }

    @Test
    public void sortCorrectlyOnCodon() {
        ReportableVariant base = ReportableVariantTestFactory.create();
        ReportableVariant variant1 = ImmutableReportableVariant.builder().from(base).canonicalHgvsCodingImpact("c.-300T>A").build();
        ReportableVariant variant2 = ImmutableReportableVariant.builder().from(base).canonicalHgvsCodingImpact("c.4000T>A").build();
        ReportableVariant variant3 = ImmutableReportableVariant.builder().from(base).canonicalHgvsCodingImpact("c.500T>A").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2, variant3);

        List<ReportableVariant> sortedVariants = SomaticVariants.sort(variants);

        assertEquals(variant1, sortedVariants.get(0));
        assertEquals(variant3, sortedVariants.get(1));
        assertEquals(variant2, sortedVariants.get(2));
    }

    @Test
    public void canExtractMSIgenes() {
        ReportableVariant base = ReportableVariantTestFactory.create();
        ReportableVariant variant1 = ImmutableReportableVariant.builder().from(base).gene("MLH1").build();
        ReportableVariant variant2 = ImmutableReportableVariant.builder().from(base).gene("MSH2").build();
        ReportableVariant variant3 = ImmutableReportableVariant.builder().from(base).gene("MSH6").build();
        ReportableVariant variant4 = ImmutableReportableVariant.builder().from(base).gene("PMS2").build();
        ReportableVariant variant5 = ImmutableReportableVariant.builder().from(base).gene("EPCAM").build();
        ReportableVariant variant6 = ImmutableReportableVariant.builder().from(base).gene("BRAF").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2, variant3, variant4, variant5, variant6);

        assertEquals(5, SomaticVariants.determineMSIgenes(variants).size());
    }

    @Test
    public void canExtractHRDgenes() {
        ReportableVariant base = ReportableVariantTestFactory.create();
        ReportableVariant variant1 = ImmutableReportableVariant.builder().from(base).gene("BRCA1").build();
        ReportableVariant variant2 = ImmutableReportableVariant.builder().from(base).gene("BRCA2").build();
        ReportableVariant variant3 = ImmutableReportableVariant.builder().from(base).gene("PALB2").build();
        ReportableVariant variant4 = ImmutableReportableVariant.builder().from(base).gene("RAD51B").build();
        ReportableVariant variant5 = ImmutableReportableVariant.builder().from(base).gene("RAD51C").build();
        ReportableVariant variant6 = ImmutableReportableVariant.builder().from(base).gene("BRAF").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2, variant3, variant4, variant5, variant6);

        assertEquals(5, SomaticVariants.determineHRDgenes(variants).size());
    }
}