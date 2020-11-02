package com.hartwig.hmftools.protect.variants;

import static com.hartwig.hmftools.protect.variants.ReportableVariantFactory.reportableGermlineVariants;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.germline.ImmutableReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantFactoryTest {

    @Test
    public void testSingleHitNotIncluded() {
        List<ReportableVariant> victims = reportableGermlineVariants(Lists.newArrayList(create("Gene", false)), Collections.emptySet());
        assertEquals(0, victims.size());
    }

    @Test
    public void testBiallelic() {
        List<ReportableVariant> victims = reportableGermlineVariants(Lists.newArrayList(create("Gene", true)), Collections.emptySet());
        assertEquals(1, victims.size());
        assertEquals("Gene", victims.get(0).gene());
    }

    @Test
    public void testKit() {
        List<ReportableVariant> victims = reportableGermlineVariants(Lists.newArrayList(create("KIT", false)), Collections.emptySet());
        assertEquals(1, victims.size());
        assertEquals("KIT", victims.get(0).gene());
    }

    @Test
    public void testDoubleGermlineHit() {
        List<ReportableVariant> victims =
                reportableGermlineVariants(Lists.newArrayList(create(1, "Gene", false), create(2, "Gene", false)), Collections.emptySet());
        assertEquals(2, victims.size());
        assertEquals("Gene", victims.get(0).gene());
        assertEquals("Gene", victims.get(1).gene());
    }

    @Test
    public void testSomaticHit() {
        List<ReportableVariant> victims =
                reportableGermlineVariants(Lists.newArrayList(create("Gene", false)), Collections.singleton("Gene"));
        assertEquals(1, victims.size());
        assertEquals("Gene", victims.get(0).gene());
    }

    @NotNull
    private static ReportableGermlineVariant create(@NotNull String gene, boolean biallelic) {
        return create(1, gene, biallelic);
    }

    @NotNull
    private static ReportableGermlineVariant create(int position, @NotNull String gene, boolean biallelic) {
        return ImmutableReportableGermlineVariant.builder()
                .gene(gene)
                .biallelic(biallelic)
                .chromosome("1")
                .position(position)
                .ref("C")
                .alt("G")
                .codingEffect(CodingEffect.MISSENSE)
                .hgvsCoding("coding")
                .hgvsProtein("protein")
                .alleleReadCount(1)
                .totalReadCount(10)
                .adjustedVaf(0.4)
                .adjustedCopyNumber(2)
                .build();
    }
}
