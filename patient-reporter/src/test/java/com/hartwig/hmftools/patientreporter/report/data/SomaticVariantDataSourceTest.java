package com.hartwig.hmftools.patientreporter.report.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SomaticVariantDataSourceTest {

    @Test
    public void canExtractCodingFromHGVSCodingImpactField() {
        assertEquals(927, SomaticVariantDataSource.extractCodonField("c.927+1G>A"));
        assertEquals(1799, SomaticVariantDataSource.extractCodonField("c.1799T>A"));
        assertEquals(423, SomaticVariantDataSource.extractCodonField("c.423_427delCCCTG"));
        assertEquals(8390, SomaticVariantDataSource.extractCodonField("c.8390delA"));
    }

    @Test
    public void sortCorrectlyOnCodon() {
        ReportableVariant variant1 = testBuilder().hgvsCodingImpact("c.300T>A").build();
        ReportableVariant variant2 = testBuilder().hgvsCodingImpact("c.4000T>A").build();
        ReportableVariant variant3 = testBuilder().hgvsCodingImpact("c.500T>A").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2, variant3);

        List<ReportableVariant> sortedVariants = SomaticVariantDataSource.sort(variants);

        assertEquals(variant1, sortedVariants.get(0));
        assertEquals(variant3, sortedVariants.get(1));
        assertEquals(variant2, sortedVariants.get(2));
    }

    @NotNull
    private static ImmutableReportableVariant.Builder testBuilder() {

        return ImmutableReportableVariant.builder()
                .gene("XXX")
                .isDrupActionable(true)
                .hgvsCodingImpact("c.1A>T")
                .hgvsProteinImpact("p.V1G")
                .hotspot(Hotspot.HOTSPOT)
                .clonality(Clonality.CLONAL)
                .alleleReadCount(1)
                .totalReadCount(1)
                .adjustedCopyNumber(1)
                .minorAllelePloidy(1)
                .biallelic(true)
                .adjustedVAF(0D)
                .driverCategory(DriverCategory.ONCO)
                .driverLikelihood(1D)
                .SomaticOrGermline("somatic")
                .notifyClinicalGeneticus(false);
    }
}