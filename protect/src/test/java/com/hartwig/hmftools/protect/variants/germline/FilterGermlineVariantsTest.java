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
    public void checkForGermlinesGenesReportedBialleic() {
        GermlineReportingModel germlineReportingModelBiallic = ProtectTestFactory.createTestGermlineGenesReportingBialleic();

        List<ReportableGermlineVariant> germlineVariantsNotPresentInTumor = createTestGermlineVariantsNotPresentInTumor();
        List<DriverSomaticVariant> driverSomaticVariants = Lists.newArrayList();
        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList();
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> reportableGeneDisruptions = Lists.newArrayList();
        List<DriverGermlineVariant> filteredGermlineVariantMatchNotPresentInTumor =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineReportingModelBiallic,
                        germlineVariantsNotPresentInTumor,
                        driverSomaticVariants,
                        reportableGainLosses,
                        reportableHomozygousDisruptions,
                        reportableGeneDisruptions);
        assertEquals(0, filteredGermlineVariantMatchNotPresentInTumor.size());

        List<ReportableGermlineVariant> germlineVariantsMultiplePresentInTumor = createTestMultipleGermlineVariantsPresentInTumor();
        List<DriverGermlineVariant> filteredGermlineVariantMatchMultipleGermline =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineReportingModelBiallic,
                        germlineVariantsMultiplePresentInTumor,
                        driverSomaticVariants,
                        reportableGainLosses,
                        reportableHomozygousDisruptions,
                        reportableGeneDisruptions);
        assertEquals(2, filteredGermlineVariantMatchMultipleGermline.size());

        List<ReportableGermlineVariant> germlineVariantsPresentInTumor = createTestGermlineVariantsPresentInTumor();
        List<DriverSomaticVariant> driverWithSomaticVariants = createSomaticVariantListForGene(ProtectTestFactory.ONCOGENE);
        List<DriverGermlineVariant> filteredGermlineVariantMatchSomatic = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModelBiallic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
        assertEquals(1, filteredGermlineVariantMatchSomatic.size());

        List<ReportableGainLoss> Reportableloss = createReportableLoss();
        List<DriverGermlineVariant> filteredGermlineVariantMatchLoss = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModelBiallic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                Reportableloss,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
        assertEquals(1, filteredGermlineVariantMatchLoss.size());

        List<ReportableGainLoss> reportableGain = createReportableGain();
        List<DriverGermlineVariant> filteredGermlineVariantMatchGain = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModelBiallic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                reportableGain,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
        assertEquals(1, filteredGermlineVariantMatchGain.size());

        List<ReportableHomozygousDisruption> homDisruption = createReportableHomozygousDisruptionListForGene(ProtectTestFactory.ONCOGENE);
        List<DriverGermlineVariant> filteredGermlineVariantMatchHomDisruption = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModelBiallic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                Reportableloss,
                homDisruption,
                reportableGeneDisruptions);
        assertEquals(1, filteredGermlineVariantMatchHomDisruption.size());

        List<ReportableGeneDisruption> disruption = createReportableGeneDisruptionListForGene(ProtectTestFactory.ONCOGENE);
        List<DriverGermlineVariant> filteredGermlineVariantMatchDisruption = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModelBiallic,
                germlineVariantsPresentInTumor,
                driverWithSomaticVariants,
                Reportableloss,
                reportableHomozygousDisruptions,
                disruption);
        assertEquals(1, filteredGermlineVariantMatchDisruption.size());

        List<ReportableGermlineVariant> germlineVariantsBiallelic = createTestGermlineVariantBiallelic(true);
        List<DriverGermlineVariant> filteredGermlineVariantMatchBialleic = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineReportingModelBiallic,
                germlineVariantsBiallelic,
                driverWithSomaticVariants,
                reportableGainLosses,
                reportableHomozygousDisruptions,
                reportableGeneDisruptions);
        assertEquals(1, filteredGermlineVariantMatchBialleic.size());

    }

    @Test
    public void checkForGermlinesGenesReportedMonoallelic() {

        GermlineReportingModel germlineReportingModelMonoallelic = ProtectTestFactory.createTestGermlineGenesReportingMonoalleic();

        List<ReportableGermlineVariant> germlineVariantsPresentInTumor = createTestGermlineVariantsPresentInTumor();
        List<DriverSomaticVariant> driverSomaticVariants = Lists.newArrayList();
        List<ReportableGainLoss> reportableGainLosses = Lists.newArrayList();
        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> reportableGeneDisruptions = Lists.newArrayList();
        List<DriverGermlineVariant> filteredGermlineVariantMatchNotPresentInTumor =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineReportingModelMonoallelic,
                        germlineVariantsPresentInTumor,
                        driverSomaticVariants,
                        reportableGainLosses,
                        reportableHomozygousDisruptions,
                        reportableGeneDisruptions);
        assertEquals(1, filteredGermlineVariantMatchNotPresentInTumor.size());

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
}