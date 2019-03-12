package com.hartwig.hmftools.patientreporter.report.data;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class SomaticVariantDataSourceTest {

    @Test
    public void canExtractCodingFromHGVSCodingImpactField() {
        assertEquals(927, SomaticVariantDataSource.extractCodonField("c.927+1G>A"));
        assertEquals(1799, SomaticVariantDataSource.extractCodonField("c.1799T>A"));
        assertEquals(423, SomaticVariantDataSource.extractCodonField("c.423_427delCCCTG"));
        assertEquals(8390, SomaticVariantDataSource.extractCodonField("c.8390delA"));
    }
}