package com.hartwig.hmftools.patientreporter.variants.germline;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ImmutableChordAnalysis;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FilterGermlineVariantsTest {

    private static final String ONCOGENE = "ONCO";
    private static final String TSG = "TSG";

    private static final DriverGeneView TEST_DRIVER_GENE_VIEW = PatientReporterTestFactory.createTestDriverGeneView(ONCOGENE, TSG);

    @NotNull
    private static ChordAnalysis createChordAnalysisNoHRD() {
        double brca1Value = 0.3;
        double brca2Value = 0.19;

        return ImmutableChordAnalysis.builder()
                .noneValue(1 - (brca1Value + brca2Value))
                .BRCA1Value(brca1Value)
                .BRCA2Value(brca2Value)
                .hrdValue(brca1Value + brca2Value)
                .predictedResponseValue(brca1Value + brca2Value > 0.5 ? true : false)
                .build();
    }

    @NotNull
    private static ChordAnalysis createChordAnalysisHRD() {
        double brca1Value = 0.6;
        double brca2Value = 0.2;

        return ImmutableChordAnalysis.builder()
                .noneValue(1 - (brca1Value + brca2Value))
                .BRCA1Value(brca1Value)
                .BRCA2Value(brca2Value)
                .hrdValue(brca1Value + brca2Value)
                .predictedResponseValue(brca1Value + brca2Value > 0.5 ? true : false)
                .build();
    }

    @Test
    public void checkForGermlineGenesReportedONCO() {
        GermlineReportingModel germlineReportingModel = PatientReporterTestFactory.createTestGermlineGenesReporting();

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        List<SomaticVariant> somaticVariants = Lists.newArrayList();

        List<GermlineVariant> germlineVariants = createTestGermlineVariantsONCOGene();
        List<InterpretGermlineVariant> filteredGermlineVariantMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariants,
                TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbers,
                somaticVariants, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantMatch.size());
    }

    @Test
    public void checkForGermlineGenesReportedTSG() {
        GermlineReportingModel germlineReportingModel = PatientReporterTestFactory.createTestGermlineGenesReporting();

        List<GermlineVariant> germlineVariantsMatch = createTestGermlineVariantsTSGGene(true, 1);
        List<GeneCopyNumber> geneCopyNumbersMatch = createCopyNumberListForTSG(1);
        List<SomaticVariant> variantsMatch = createEnrichedListForGene(TSG);
        List<InterpretGermlineVariant> filteredGermlineVariantMatch =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariantsMatch, TEST_DRIVER_GENE_VIEW,
                        germlineReportingModel,
                        geneCopyNumbersMatch,
                        variantsMatch, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantMatch.size()); // all three options matched

        List<GermlineVariant> germlineVariantsNonMatchBiallelic = createTestGermlineVariantsTSGGene(false, 1);
        List<GeneCopyNumber> geneCopyNumbersNonMatchBiallelic = createCopyNumberListForTSG(1);
        List<SomaticVariant> variantsNonMatchBiallelic = createEnrichedListForGene(TSG);
        List<InterpretGermlineVariant> filteredGermlineVariantNonMatchBiallelic = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatchBiallelic, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatchBiallelic,
                variantsNonMatchBiallelic, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantNonMatchBiallelic.size()); // match copy number and variant

        List<GermlineVariant> germlineVariantsNonMatchVariant = createTestGermlineVariantsTSGGene(true, 1);
        List<GeneCopyNumber> geneCopyNumbersNonMatchVariant = createCopyNumberListForTSG(1);
        List<SomaticVariant> variantsNonMatchVariant = createEnrichedListForGene("AAAA");
        List<InterpretGermlineVariant> filteredGermlineVariantNonMatchVariant = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatchVariant, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatchVariant,
                variantsNonMatchVariant, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantNonMatchVariant.size()); // match biallelic and copy number

        List<GermlineVariant> germlineVariantsNonMatchCopy = createTestGermlineVariantsTSGGene(true, 1);
        List<GeneCopyNumber> geneCopyNumbersNonMatchCopy = createCopyNumberListForTSG(2);
        List<SomaticVariant> variantsNonMatchCopy = createEnrichedListForGene(TSG);
        List<InterpretGermlineVariant> filteredGermlineVariantNonMatchCopy = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatchCopy, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatchCopy,
                variantsNonMatchCopy, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantNonMatchCopy.size()); // match biallelic and variant

        List<GermlineVariant> germlineVariantsNonMatch = createTestGermlineVariantsTSGGene(false, 1);
        List<GeneCopyNumber> geneCopyNumbersNonMatch = createCopyNumberListForTSG(2);
        List<SomaticVariant> variantsNonMatch = createEnrichedListForGene("AAAA");
        List<InterpretGermlineVariant> filteredGermlineVariantNonMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatch, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersNonMatch,
                variantsNonMatch, createChordAnalysisNoHRD());
        assertEquals(0, filteredGermlineVariantNonMatch.size()); // all option failed

        List<GermlineVariant> germlineVariantsOptionBiallelic = createTestGermlineVariantsTSGGene(true, 1);
        List<GeneCopyNumber> geneCopyNumbersOptionBiallelic = createCopyNumberListForTSG(2);
        List<SomaticVariant> variantsOptionBiallelic = createEnrichedListForGene("AAAA");
        List<InterpretGermlineVariant> filteredGermlineVariantOptionBiallelic = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionBiallelic, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersOptionBiallelic,
                variantsOptionBiallelic, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantOptionBiallelic.size()); // only match biallelic

        List<GermlineVariant> germlineVariantsOptionVariant = createTestGermlineVariantsTSGGene(false, 1);
        List<GeneCopyNumber> geneCopyNumbersOptionVariant = createCopyNumberListForTSG(2);
        List<SomaticVariant> variantsOptionVariant = createEnrichedListForGene(TSG);
        List<InterpretGermlineVariant> filteredGermlineVariantOptionVariant = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionVariant, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersOptionVariant,
                variantsOptionVariant, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantOptionVariant.size()); // only match variant

        List<GermlineVariant> germlineVariantsOptionCopyNumberPartialLoss = createTestGermlineVariantsTSGGene(false, 1);
        List<GeneCopyNumber> geneCopyNumbersCopyNumberPartialLoss = createCopyNumberListForTSG(1);
        List<SomaticVariant> variantsOptionCopyNumberPartialLoss = createEnrichedListForGene("AAAA");
        List<InterpretGermlineVariant> filteredGermlineVariantOptionCopyNumberPartialLoss = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionCopyNumberPartialLoss, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersCopyNumberPartialLoss,
                variantsOptionCopyNumberPartialLoss, createChordAnalysisNoHRD());
        assertEquals(0, filteredGermlineVariantOptionCopyNumberPartialLoss.size()); // only match copy number

        List<GermlineVariant> germlineVariantsOptionCopyNumberFullLoss = createTestGermlineVariantsTSGGene(false, 2);
        List<GeneCopyNumber> geneCopyNumbersCopyNumberFullLoss = createCopyNumberListForTSG(1);
        List<SomaticVariant> variantsOptionCopyNumberFullLoss = createEnrichedListForGene("AAAA");
        List<InterpretGermlineVariant> filteredGermlineVariantOptionCopyNumberFullLoss = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionCopyNumberFullLoss, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersCopyNumberFullLoss,
                variantsOptionCopyNumberFullLoss, createChordAnalysisNoHRD());
        assertEquals(1, filteredGermlineVariantOptionCopyNumberFullLoss.size()); // only match copy number

        List<GermlineVariant> germlineVariantsOptionCopyNumberHRD = createTestGermlineVariantsTSGGene(false, 2);
        List<GeneCopyNumber> geneCopyNumbersCopyNumberHRD = createCopyNumberListForTSG(4);
        List<SomaticVariant> variantsOptionCopyNumberHRD = createEnrichedListForGene("AAAA");
        List<InterpretGermlineVariant> filteredGermlineVariantOptionCopyNumberHRD = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionCopyNumberHRD, TEST_DRIVER_GENE_VIEW,
                germlineReportingModel,
                geneCopyNumbersCopyNumberHRD,
                variantsOptionCopyNumberHRD, createChordAnalysisHRD());
        assertEquals(1, filteredGermlineVariantOptionCopyNumberHRD.size()); // only match HRD
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsONCOGene() {
        return Lists.newArrayList(PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(ONCOGENE).build());
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsTSGGene(boolean biallelicFilter, double copyNumber) {
        return Lists.newArrayList(PatientReporterTestFactory.createTestGermlineVariantBuilder()
                .gene(TSG)
                .biallelic(biallelicFilter)
                .adjustedCopyNumber(copyNumber)
                .build());
    }

    @NotNull
    private static List<GeneCopyNumber> createCopyNumberListForTSG(int minCopyNumber) {
        return Lists.newArrayList(PatientReporterTestFactory.createTestCopyNumberBuilder().gene(TSG).minCopyNumber(minCopyNumber).build());
    }

    @NotNull
    private static List<SomaticVariant> createEnrichedListForGene(@NotNull String gene) {
        return Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(gene).build());
    }
}