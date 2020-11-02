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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FilterGermlineVariantsTest {

    @Test
    public void checkForGermlineGenesReportedBiallelic() {
        GermlineReportingModel germlineReportingModelBiallelic =
                ProtectTestFactory.createTestGermlineGenesReporting(true, false, true, true, null, null);

        List<ReportableGermlineVariant> germlineVariantsNotPresentInTumor = createTestGermlineVariantsNotPresentInTumor();
        List<DriverSomaticVariant> driverSomaticVariants = Lists.newArrayList();
        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList();
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> reportableGeneDisruptions = Lists.newArrayList();
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsNotPresentInTumor,
                driverSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);

        List<ReportableGermlineVariant> germlineVariantsMultiplePresentInTumor = createTestMultipleGermlineVariantsPresentInTumor();
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsMultiplePresentInTumor,
                driverSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);

        List<ReportableGermlineVariant> germlineVariantsPresentInTumor = createTestGermlineVariantsPresentInTumor();
        List<DriverSomaticVariant> driverWithSomaticVariants = createSomaticVariantListForGene(ProtectTestFactory.ONCOGENE);
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);

        List<ReportableGainLoss> Reportableloss = createReportableLoss();
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                Reportableloss,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);

        List<ReportableGainLoss> reportableGain = createReportableGain();
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                reportableGain,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);

        List<ReportableHomozygousDisruption> homDisruption = createReportableHomozygousDisruptionListForGene(ProtectTestFactory.ONCOGENE);
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                Reportableloss,
                homDisruption,
                reportableGeneDisruptions);

        List<ReportableGeneDisruption> disruption = createReportableGeneDisruptionListForGene(ProtectTestFactory.ONCOGENE);
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                Reportableloss,
                reportableHomozygousDisruptions,
                disruption);

        List<ReportableGermlineVariant> germlineVariantsBiallelic = createTestGermlineVariantBiallelic(true);
        assertFilter(0,
                germlineReportingModelBiallelic,
                germlineVariantsBiallelic,
                driverWithSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);

    }

    private static void assertFilter(int sizeExpected, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull List<ReportableGermlineVariant> reportableGermlineVariants, @NotNull List<DriverSomaticVariant> driverSomaticVariant,
            @NotNull List<ReportableGainLoss> reportableGainLosses,
            @NotNull List<ReportableHomozygousDisruption> reportableHomozygousDisruptions,
            @NotNull List<ReportableGeneDisruption> reportableGeneDisruptions) {

        List<DriverGermlineVariant> filteredGermlineVariants = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModel,
                reportableGermlineVariants,
                driverSomaticVariant,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
        assertEquals(sizeExpected, filteredGermlineVariants.size());
    }

    @Test
    public void checkForGermlineGenesReportedMonoallelic() {
        GermlineReportingModel germlineReportingModelMonoallelic =
                ProtectTestFactory.createTestGermlineGenesReporting(true, false, false, false, "ABC", null);

        List<ReportableGermlineVariant> germlineVariantsPresentInTumor = createTestGermlineVariantsPresentInTumor();
        List<DriverSomaticVariant> driverSomaticVariants = Lists.newArrayList();
        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList();
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> reportableGeneDisruptions = Lists.newArrayList();
        assertFilter(1, germlineReportingModelMonoallelic, germlineVariantsPresentInTumor,
                driverSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
    }

    @NotNull
    private static List<ReportableGermlineVariant> createTestGermlineVariantsNotPresentInTumor() {
        return Lists.newArrayList(ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .adjustedVaf(0.1)
                .adjustedCopyNumber(1D)
                .hgvsProtein(Strings.EMPTY)
                .build());
    }

    @NotNull
    private static List<ReportableGermlineVariant> createTestMultipleGermlineVariantsPresentInTumor() {
        ReportableGermlineVariant reportableGermlineVariant1 = ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .adjustedVaf(0.6)
                .adjustedCopyNumber(1D)
                .hgvsProtein("ABC")
                .build();

        ReportableGermlineVariant reportableGermlineVariant2 = ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .adjustedVaf(0.6)
                .adjustedCopyNumber(1D)
                .hgvsProtein("DEF")
                .build();

        return Lists.newArrayList(reportableGermlineVariant1, reportableGermlineVariant2);
    }

    @NotNull
    private static ReportableGermlineVariant createTestGermlineVariants() {
        return ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .adjustedVaf(0.6)
                .adjustedCopyNumber(1D)
                .hgvsProtein("ABC")
                .build();
    }

    @NotNull
    private static List<ReportableGermlineVariant> createTestGermlineVariantsPresentInTumor() {
        return Lists.newArrayList(ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.ONCOGENE)
                .adjustedVaf(0.6)
                .adjustedCopyNumber(1D)
                .hgvsProtein("ABC")
                .build());
    }

    @NotNull
    private static List<DriverSomaticVariant> createSomaticVariantListForGene(@NotNull String gene) {
        return Lists.newArrayList(ImmutableDriverSomaticVariant.builder()
                .variant(ProtectTestFactory.createTestSomaticVariantBuilder().gene(gene).build())
                .driverLikelihood(1D)
                .build());
    }

    @NotNull
    private static List<ReportableGermlineVariant> createTestGermlineVariantBiallelic(boolean biallelic) {
        return Lists.newArrayList(ProtectTestFactory.createTestGermlineVariantBuilder()
                .gene(ProtectTestFactory.TSG)
                .biallelic(biallelic)
                .adjustedCopyNumber(1D)
                .adjustedVaf(biallelic ? 1D : 0.5)
                .hgvsProtein(Strings.EMPTY)
                .build());
    }

    @NotNull
    private static List<ReportableGainLoss> createReportableLoss() {
        return Lists.newArrayList(ProtectTestFactory.createTestReportableGainLossBuilder()
                .gene(ProtectTestFactory.TSG)
                .copies(0)
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .build());
    }

    @NotNull
    private static List<ReportableGainLoss> createReportableGain() {
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

    @Test
    public void checkForGermlineGenesReportedMonoallelicSpecificVariant() {
        GermlineReportingModel germlineReportingModelMonoallelic =
                ProtectTestFactory.createTestGermlineGenesReporting(true, false, false, false, "ABC", null);

        List<ReportableGermlineVariant> germlineVariantsPresentInTumor = createTestGermlineVariantsPresentInTumor();
        List<DriverSomaticVariant> driverSomaticVariants = Lists.newArrayList();
        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList();
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> reportableGeneDisruptions = Lists.newArrayList();

        assertFilter(1, germlineReportingModelMonoallelic,
                germlineVariantsPresentInTumor,
                driverSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
    }
}