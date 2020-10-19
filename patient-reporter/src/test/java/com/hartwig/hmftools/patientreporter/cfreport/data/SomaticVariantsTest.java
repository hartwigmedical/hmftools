package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.protect.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
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
        ReportableVariant variant1 = createTestReportableVariantBuilder().canonicalHgvsCodingImpact("c.-300T>A").build();
        ReportableVariant variant2 = createTestReportableVariantBuilder().canonicalHgvsCodingImpact("c.4000T>A").build();
        ReportableVariant variant3 = createTestReportableVariantBuilder().canonicalHgvsCodingImpact("c.500T>A").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2, variant3);

        List<ReportableVariant> sortedVariants = SomaticVariants.sort(variants);

        assertEquals(variant1, sortedVariants.get(0));
        assertEquals(variant3, sortedVariants.get(1));
        assertEquals(variant2, sortedVariants.get(2));
    }

    @NotNull
    private static ImmutableReportableVariant.Builder createTestReportableVariantBuilder() {
        return ImmutableReportableVariant.builder()
                .gene(Strings.EMPTY)
                .position(0)
                .chromosome(Strings.EMPTY)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .gDNA(Strings.EMPTY)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .alleleReadCount(0)
                .totalReadCount(0)
                .alleleCopyNumber(0D)
                .totalCopyNumber(0)
                .biallelic(false)
                .driverLikelihood(0D)
                .notifyClinicalGeneticist(false);
    }
}