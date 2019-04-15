package com.hartwig.hmftools.patientreporter.report.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.germline.ImmutableGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.ImmutableReportableSomaticVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableSomaticVariant;

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
        ReportableSomaticVariant variant1 = testBuilder().hgvsCodingImpact("c.300T>A").build();
        ReportableSomaticVariant variant2 = testBuilder().hgvsCodingImpact("c.4000T>A").build();
        ReportableSomaticVariant variant3 = testBuilder().hgvsCodingImpact("c.500T>A").build();

        List<ReportableSomaticVariant> variants = Lists.newArrayList(variant1, variant2, variant3);

        List<ReportableSomaticVariant> sortedVariants = SomaticVariantDataSource.sort(variants);

        assertEquals(variant1, sortedVariants.get(0));
        assertEquals(variant3, sortedVariants.get(1));
        assertEquals(variant2, sortedVariants.get(2));
    }

    @NotNull
    private static ImmutableReportableSomaticVariant.Builder testBuilder() {
        List<GermlineVariant> germlineVariant = Lists.newArrayList(ImmutableGermlineVariant.builder()
                .passFilter(true)
                .gene("BRCA2")
                .hgvsCodingImpact("c.5946delT")
                .hgvsProteinImpact("p.Ser1982fs")
                .totalReadCount(112)
                .alleleReadCount(67)
                .germlineStatus("HET")
                .adjustedCopyNumber(3D)
                .adjustedVAF(1.0)
                .minorAllelePloidy(1D)
                .biallelic(true)
                .build());

        return ImmutableReportableSomaticVariant.builder()
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
                .germlineVariant(germlineVariant)
                .SomaticOrGermline("somatic");
    }
}