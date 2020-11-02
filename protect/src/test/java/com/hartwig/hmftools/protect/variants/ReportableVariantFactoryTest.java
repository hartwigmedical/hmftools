package com.hartwig.hmftools.protect.variants;

import static com.hartwig.hmftools.protect.variants.ReportableVariantFactory.reportableGermlineVariants;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.germline.ImmutableReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.protect.ProtectTestFactory;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantFactoryTest {

    private static final String GENE = "Gene";

    @Test
    public void singleHitNotIncluded() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineModel(GENE, true, null);
        List<ReportableGermlineVariant> variants = Lists.newArrayList(create(GENE, false));

        List<ReportableVariant> victims = reportableGermlineVariants(variants, Collections.emptySet(), germlineReportingModel);
        assertEquals(0, victims.size());
    }

    @Test
    public void singleHitIncludedWhenAllowed() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineModel(GENE, false, null);
        List<ReportableGermlineVariant> variants = Lists.newArrayList(create(GENE, false));

        List<ReportableVariant> victims = reportableGermlineVariants(variants, Collections.emptySet(), germlineReportingModel);
        assertEquals(1, victims.size());
        assertEquals(GENE, victims.get(0).gene());
    }

    @Test
    public void biallelicIncluded() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineModel(GENE, true, null);
        List<ReportableGermlineVariant> variants = Lists.newArrayList(create(GENE, true));

        List<ReportableVariant> victims = reportableGermlineVariants(variants, Collections.emptySet(), germlineReportingModel);
        assertEquals(1, victims.size());
        assertEquals(GENE, victims.get(0).gene());
    }

    @Test
    public void biallelicExcludedWhenNotInGermlineList() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createEmptyGermlineReportingModel();
        List<ReportableGermlineVariant> variants = Lists.newArrayList(create(GENE, true));

        List<ReportableVariant> victims = reportableGermlineVariants(variants, Collections.emptySet(), germlineReportingModel);
        assertEquals(0, victims.size());
    }

    @Test
    public void biallelicExcludedWhenNotInTumor() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineModel(GENE, true, null);
        List<ReportableGermlineVariant> variants = Lists.newArrayList(create(GENE, 1, "protein", 0.1, true));

        List<ReportableVariant> victims = reportableGermlineVariants(variants, Collections.emptySet(), germlineReportingModel);
        assertEquals(0, victims.size());
    }

    @Test
    public void singleHitIncludedOnExtraGermlineHit() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineModel(GENE, true, null);
        List<ReportableGermlineVariant> variants = Lists.newArrayList(create(GENE, 1, false), create(GENE, 2, false));

        List<ReportableVariant> victims = reportableGermlineVariants(variants, Collections.emptySet(), germlineReportingModel);
        assertEquals(2, victims.size());
        assertEquals(GENE, victims.get(0).gene());
        assertEquals(GENE, victims.get(1).gene());
    }

    @Test
    public void singleHitIncludedOnExtraSomaticHit() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineModel(GENE, true, null);
        List<ReportableGermlineVariant> variants = Lists.newArrayList(create(GENE, false));

        List<ReportableVariant> victims = reportableGermlineVariants(variants, Collections.singleton(GENE), germlineReportingModel);
        assertEquals(1, victims.size());
        assertEquals(GENE, victims.get(0).gene());
    }

    @Test
    public void exclusiveHgvsProteinFilterWorks() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineModel(GENE, true, "proteinMatch");
        List<ReportableGermlineVariant> variantMatch = Lists.newArrayList(create(GENE, 1, "proteinMatch", 0.4, false));

        List<ReportableVariant> victimsMatch = reportableGermlineVariants(variantMatch, Collections.emptySet(), germlineReportingModel);
        assertEquals(1, victimsMatch.size());
        assertEquals(GENE, victimsMatch.get(0).gene());

        List<ReportableGermlineVariant> variantsNonMatch = Lists.newArrayList(create(GENE, 1, "weirdProtein", 0.4, false));
        List<ReportableVariant> victimsNonMatch =
                reportableGermlineVariants(variantsNonMatch, Collections.emptySet(), germlineReportingModel);
        assertEquals(0, victimsNonMatch.size());
    }

    @NotNull
    private static ReportableGermlineVariant create(@NotNull String gene, boolean biallelic) {
        return create(gene, 1, biallelic);
    }

    @NotNull
    private static ReportableGermlineVariant create(@NotNull String gene, int position, boolean biallelic) {
        return create(gene, position, "protein", 0.4, biallelic);
    }

    @NotNull
    private static ReportableGermlineVariant create(@NotNull String gene, int position, @NotNull String hgvsProtein, double adjustedVaf,
            boolean biallelic) {
        return ImmutableReportableGermlineVariant.builder()
                .gene(gene)
                .chromosome("1")
                .biallelic(biallelic)
                .position(position)
                .ref("C")
                .alt("G")
                .codingEffect(CodingEffect.MISSENSE)
                .hgvsCoding("coding")
                .hgvsProtein(hgvsProtein)
                .alleleReadCount(1)
                .totalReadCount(10)
                .adjustedVaf(adjustedVaf)
                .adjustedCopyNumber(2)
                .build();
    }
}
