package com.hartwig.hmftools.protect.variants.germline;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.protect.ProtectTestFactory;
import com.hartwig.hmftools.protect.homozygousdisruption.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.structural.ReportableGeneDisruption;
import com.hartwig.hmftools.protect.variants.somatic.DriverSomaticVariant;
import com.hartwig.hmftools.protect.variants.somatic.ImmutableDriverSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FilterGermlineVariantsTest {

    @Test
    public void checkForGermlineGenesReportedONCO() {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineGenesReporting();

        List<DriverSomaticVariant> driverSomaticVariants = Lists.newArrayList();
        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList();
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> reportableGeneDisruptions = Lists.newArrayList();

        List<ReportableGermlineVariant> germlineVariantsNotPresentInTumor = createTestGermlineVariantsOncogeneNotPresentInTumor();
        List<DriverGermlineVariant> filteredGermlineVariantMatchNotPresentInTumor =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineReportingModel,
                        germlineVariantsNotPresentInTumor,
                        driverSomaticVariants,
                        reportableGainLosses,
                        reportableHomozygousDisruptions,
                        reportableGeneDisruptions);
        assertEquals(0, filteredGermlineVariantMatchNotPresentInTumor.size());

        List<ReportableGermlineVariant> germlineVariantsPresentInTumor = createTestGermlineVariantsOncogenePresentInTumor();
        List<DriverGermlineVariant> filteredGermlineVariantMatchPresentInTumor = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModel,
                germlineVariantsPresentInTumor,
                driverSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
        assertEquals(1, filteredGermlineVariantMatchPresentInTumor.size());
    }

    @Test
    public void checkForGermlineGenesReportedTSG() {
        List<ReportableGermlineVariant> germlineVariantsBiallelic = createTestGermlineVariantListTSG(true);
        List<ReportableGermlineVariant> germlineVariantsNonBiallelic = createTestGermlineVariantListTSG(false);
        List<ReportableGermlineVariant> germlineVariantsDoubleHit = Lists.newArrayList();
        germlineVariantsDoubleHit.addAll(createTestGermlineVariantListTSG(false));
        germlineVariantsDoubleHit.addAll(createTestGermlineVariantListTSG(false));

        List<DriverSomaticVariant> variantOnTSG = createSomaticVariantListForGene(ProtectTestFactory.TSG);
        List<DriverSomaticVariant> variantOnOncogene = createSomaticVariantListForGene(ProtectTestFactory.ONCOGENE);

        List<ReportableGainLoss> lossOnTSG = createReportableLossOnTSG();
        List<ReportableGainLoss> gainOnOncogene = createReportableGainOnOncogene();

        List<ReportableHomozygousDisruption> homDisruptionOnTSG = createReportableHomozygousDisruptionListForGene(ProtectTestFactory.TSG);
        List<ReportableHomozygousDisruption> homDisruptionOnOncogene =
                createReportableHomozygousDisruptionListForGene(ProtectTestFactory.ONCOGENE);

        List<ReportableGeneDisruption> disruptionOnTSG = createReportableGeneDisruptionListForGene(ProtectTestFactory.TSG);
        List<ReportableGeneDisruption> disruptionOnOncogene = createReportableGeneDisruptionListForGene(ProtectTestFactory.ONCOGENE);

        assertFilter(true, germlineVariantsBiallelic, variantOnTSG, lossOnTSG, homDisruptionOnTSG, disruptionOnTSG);
        assertFilter(true, germlineVariantsBiallelic, variantOnOncogene, gainOnOncogene, homDisruptionOnOncogene, disruptionOnOncogene);
        assertFilter(true, germlineVariantsNonBiallelic, variantOnTSG, gainOnOncogene, homDisruptionOnOncogene, disruptionOnOncogene);
        assertFilter(true, germlineVariantsNonBiallelic, variantOnOncogene, lossOnTSG, homDisruptionOnOncogene, disruptionOnOncogene);
        assertFilter(true, germlineVariantsNonBiallelic, variantOnOncogene, gainOnOncogene, homDisruptionOnTSG, disruptionOnOncogene);
        assertFilter(true, germlineVariantsNonBiallelic, variantOnOncogene, gainOnOncogene, homDisruptionOnOncogene, disruptionOnTSG);
        assertFilter(true, germlineVariantsDoubleHit, variantOnOncogene, gainOnOncogene, homDisruptionOnOncogene, disruptionOnOncogene);
        assertFilter(false, germlineVariantsNonBiallelic, variantOnOncogene, gainOnOncogene, homDisruptionOnOncogene, disruptionOnOncogene);
    }

    private static void assertFilter(boolean expectedToPass, @NotNull List<ReportableGermlineVariant> reportableGermlineVariants,
            @NotNull List<DriverSomaticVariant> driverSomaticVariant, @NotNull List<ReportableGainLoss> reportableGainLosses,
            @NotNull List<ReportableHomozygousDisruption> reportableHomozygousDisruptions,
            @NotNull List<ReportableGeneDisruption> reportableGeneDisruptions) {
        GermlineReportingModel germlineReportingModel = ProtectTestFactory.createTestGermlineGenesReporting();

        List<DriverGermlineVariant> filteredGermlineVariants = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModel,
                reportableGermlineVariants,
                driverSomaticVariant,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
        assertEquals(expectedToPass, !filteredGermlineVariants.isEmpty());
    }

    @NotNull
    private static List<ReportableGermlineVariant> createTestGermlineVariantsOncogeneNotPresentInTumor() {
        return Lists.newArrayList(ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .adjustedVaf(0.1)
                .adjustedCopyNumber(1D)
                .build());
    }

    @NotNull
    private static List<ReportableGermlineVariant> createTestGermlineVariantsOncogenePresentInTumor() {
        return Lists.newArrayList(ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .adjustedVaf(0.6)
                .adjustedCopyNumber(1D)
                .build());
    }

    @NotNull
    private static List<ReportableGermlineVariant> createTestGermlineVariantListTSG(boolean biallelic) {
        return Lists.newArrayList(ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.TSG)
                .biallelic(biallelic)
                .adjustedCopyNumber(1D)
                .adjustedVaf(biallelic ? 1D : 0.5)
                .build());
    }

    @NotNull
    private static List<ReportableGainLoss> createReportableLossOnTSG() {
        return Lists.newArrayList(ProtectTestFactory.createTestReportableGainLossBuilder()
                .gene(ProtectTestFactory.TSG)
                .copies(0)
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .build());
    }

    @NotNull
    private static List<ReportableGainLoss> createReportableGainOnOncogene() {
        return Lists.newArrayList(ProtectTestFactory.createTestReportableGainLossBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .copies(100)
                .interpretation(CopyNumberInterpretation.GAIN)
                .build());
    }

    @NotNull
    private static List<ReportableHomozygousDisruption> createReportableHomozygousDisruptionListForGene(@NotNull String gene) {
        return Lists.newArrayList(ProtectTestFactory.createTestReportableHomozygousDisruptionBuilder().gene(gene).build());
    }

    @NotNull
    private static List<ReportableGeneDisruption> createReportableGeneDisruptionListForGene(@NotNull String gene) {
        return Lists.newArrayList(ProtectTestFactory.createTestReportableGeneDisruptionBuilder().gene(gene).build());
    }

    @NotNull
    private static List<DriverSomaticVariant> createSomaticVariantListForGene(@NotNull String gene) {
        return Lists.newArrayList(ImmutableDriverSomaticVariant.builder()
                .variant(ProtectTestFactory.createTestSomaticVariantBuilder().gene(gene).build())
                .driverLikelihood(0D)
                .build());
    }
}