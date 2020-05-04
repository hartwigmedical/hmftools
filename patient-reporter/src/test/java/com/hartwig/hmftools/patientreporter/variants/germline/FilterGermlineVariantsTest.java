package com.hartwig.hmftools.patientreporter.variants.germline;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordAnalyzer;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableGermlineVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FilterGermlineVariantsTest {

    private static final String ONCOGENE = "ONCO";
    private static final String TSG = "TSG";

    private static final DriverGeneView TEST_DRIVER_GENE_VIEW = PatientReporterTestFactory.createTestDriverGeneView(ONCOGENE, TSG);

    @Test
    public void checkForGermlineGenesReportedONCO() {
        GermlineReportingModel germlineReportingModel = PatientReporterTestFactory.createTestGermlineGenesReporting();

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        List<SomaticVariant> somaticVariants = Lists.newArrayList();

        List<GermlineVariant> germlineVariants = createTestGermlineVariantsONCOGene();
        List<ReportableGermlineVariant> filteredGermlineVariantMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariants,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbers,
                somaticVariants,
                createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantMatch.size());
    }

    @Test
    public void checkForGermlineGenesReportedTSG() {
        GermlineReportingModel germlineReportingModel = PatientReporterTestFactory.createTestGermlineGenesReporting();

        List<GermlineVariant> germlineVariantsMatch = createTestGermlineVariantsTSGGene(true, 1);
        List<GermlineVariant> germlineVariantsNonMatch = createTestGermlineVariantsTSGGene(false, 2);

        List<GeneCopyNumber> geneCopyNumbersMatch = createCopyNumberListForTSG(1);
        List<GeneCopyNumber> geneCopyNumbersNonMatch = createCopyNumberListForTSG(2);

        List<SomaticVariant> variantsMatch = createSomaticVariantListForGene(TSG);
        List<SomaticVariant> variantsNonMatch = createSomaticVariantListForGene("AAAA");

        List<ReportableGermlineVariant> filteredGermlineVariantAllMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsMatch,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersMatch,
                variantsMatch,
                createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantAllMatch.size()); // all three options matched

        List<ReportableGermlineVariant> filteredGermlineVariantNonMatchBiallelic =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariantsNonMatch,
                        TEST_DRIVER_GENE_VIEW,
                        germlineReportingModel,
                        geneCopyNumbersMatch,
                        variantsMatch,
                        createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantNonMatchBiallelic.size()); // match copy number and variant

        List<ReportableGermlineVariant> filteredGermlineVariantNonMatchVariant = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsMatch,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersMatch,
                variantsNonMatch,
                createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantNonMatchVariant.size()); // match biallelic and copy number

        List<ReportableGermlineVariant> filteredGermlineVariantNonMatchCopy = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsMatch,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatch,
                variantsMatch,
                createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantNonMatchCopy.size()); // match biallelic and variant

        List<ReportableGermlineVariant> filteredGermlineVariantNonMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatch,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatch,
                variantsNonMatch,
                createChordAnalysisNoHRD());
        assertEquals(0, filteredGermlineVariantNonMatch.size()); // all option failed

        List<ReportableGermlineVariant> filteredGermlineVariantOptionBiallelic = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsMatch,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatch,
                variantsNonMatch,
                createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantOptionBiallelic.size()); // only match biallelic

        List<ReportableGermlineVariant> filteredGermlineVariantOptionVariant = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatch,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatch,
                variantsMatch,
                createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantOptionVariant.size()); // only match variant

        List<ReportableGermlineVariant> filteredGermlineVariantOptionCopyNumberPartialLoss =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariantsNonMatch,
                        TEST_DRIVER_GENE_VIEW,
                        germlineReportingModel,
                        geneCopyNumbersMatch,
                        variantsNonMatch,
                        createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantOptionCopyNumberPartialLoss.size()); // only match copy number

        List<ReportableGermlineVariant> filteredGermlineVariantOptionCopyNumberHRD =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariantsNonMatch,
                        TEST_DRIVER_GENE_VIEW,
                        germlineReportingModel,
                        geneCopyNumbersNonMatch,
                        variantsNonMatch,
                        createChordAnalysisHRD());
        assertEquals(1, filteredGermlineVariantOptionCopyNumberHRD.size()); // only match HRD
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsONCOGene() {
        return Lists.newArrayList(PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(ONCOGENE).build());
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsTSGGene(boolean biallelic, double adjustedCopyNumber) {
        return Lists.newArrayList(PatientReporterTestFactory.createTestGermlineVariantBuilder()
                .gene(TSG)
                .biallelic(biallelic)
                .adjustedCopyNumber(adjustedCopyNumber)
                .build());
    }

    @NotNull
    private static List<GeneCopyNumber> createCopyNumberListForTSG(int minCopyNumber) {
        return Lists.newArrayList(PatientReporterTestFactory.createTestCopyNumberBuilder().gene(TSG).minCopyNumber(minCopyNumber).build());
    }

    @NotNull
    private static List<SomaticVariant> createSomaticVariantListForGene(@NotNull String gene) {
        return Lists.newArrayList(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(gene).build());
    }

    @NotNull
    private static ChordAnalyzer createChordAnalysisNoHRD() {
        double brca1Value = 0.3;
        double brca2Value = 0.19;

        ChordAnalysis chordAnalysis = ImmutableChordAnalysis.builder()
                .noneValue(1 - (brca1Value + brca2Value))
                .BRCA1Value(brca1Value)
                .BRCA2Value(brca2Value)
                .hrdValue(brca1Value + brca2Value)
                .predictedResponseValue(brca1Value + brca2Value > 0.5)
                .build();
        return ImmutableChordAnalyzer.builder().chordAnalysis(chordAnalysis).hrdStatus(ChordStatus.HRP).build();
    }

    @NotNull
    private static ChordAnalyzer createChordAnalysisHRD() {
        double brca1Value = 0.6;
        double brca2Value = 0.2;

        ChordAnalysis chordAnalysis =  ImmutableChordAnalysis.builder()
                .noneValue(1 - (brca1Value + brca2Value))
                .BRCA1Value(brca1Value)
                .BRCA2Value(brca2Value)
                .hrdValue(brca1Value + brca2Value)
                .predictedResponseValue(brca1Value + brca2Value > 0.5)
                .build();

        return ImmutableChordAnalyzer.builder().chordAnalysis(chordAnalysis).hrdStatus(ChordStatus.HRD).build();

    }
}