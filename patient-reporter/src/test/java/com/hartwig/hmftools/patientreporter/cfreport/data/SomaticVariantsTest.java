package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.HomozygousDisruptionFactoryTest;
import com.hartwig.hmftools.common.linx.ImmutableHomozygousDisruption;
import com.hartwig.hmftools.common.purple.loader.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.common.purple.loader.ImmutableGainLoss;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;

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
        ReportableVariant variant2 = ImmutableReportableVariant.builder().from(base).gene("BRAF").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2);

        GainLoss baseGainLoss = GainsAndLossesTest.testGainLoss("1", "p.12");
        GainLoss gainLoss1 =
                ImmutableGainLoss.builder().from(baseGainLoss).gene("MSH2").interpretation(CopyNumberInterpretation.FULL_LOSS).build();
        GainLoss gainLoss2 =
                ImmutableGainLoss.builder().from(baseGainLoss).gene("MSH6").interpretation(CopyNumberInterpretation.PARTIAL_LOSS).build();
        GainLoss gainLoss3 =
                ImmutableGainLoss.builder().from(baseGainLoss).gene("EPCAM").interpretation(CopyNumberInterpretation.FULL_GAIN).build();

        List<GainLoss> gainLosses = Lists.newArrayList(gainLoss1, gainLoss2, gainLoss3);

        List<HomozygousDisruption> homozygousDisruption = Lists.newArrayList(create("PMS2"));

        assertEquals(4, SomaticVariants.determineMSIgenes(variants, gainLosses, homozygousDisruption).size());
    }

    @Test
    public void canExtractHRDgenes() {
        ReportableVariant base = ReportableVariantTestFactory.create();
        ReportableVariant variant1 = ImmutableReportableVariant.builder().from(base).gene("BRCA1").build();
        ReportableVariant variant2 = ImmutableReportableVariant.builder().from(base).gene("BRAF").build();

        List<ReportableVariant> variants = Lists.newArrayList(variant1, variant2);

        GainLoss baseGainLoss = GainsAndLossesTest.testGainLoss("1", "p.12");
        GainLoss gainLoss1 =
                ImmutableGainLoss.builder().from(baseGainLoss).gene("BRCA2").interpretation(CopyNumberInterpretation.FULL_LOSS).build();
        GainLoss gainLoss2 =
                ImmutableGainLoss.builder().from(baseGainLoss).gene("PALB2").interpretation(CopyNumberInterpretation.PARTIAL_LOSS).build();
        GainLoss gainLoss3 =
                ImmutableGainLoss.builder().from(baseGainLoss).gene("RAD51B").interpretation(CopyNumberInterpretation.FULL_GAIN).build();

        List<GainLoss> gainLosses = Lists.newArrayList(gainLoss1, gainLoss2, gainLoss3);

        List<HomozygousDisruption> homozygousDisruption = Lists.newArrayList(create("RAD51C"));


        assertEquals(4, SomaticVariants.determineHRDgenes(variants, gainLosses, homozygousDisruption).size());
    }

    @NotNull
    private static HomozygousDisruption create(@NotNull String gene) {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .build();
    }
}