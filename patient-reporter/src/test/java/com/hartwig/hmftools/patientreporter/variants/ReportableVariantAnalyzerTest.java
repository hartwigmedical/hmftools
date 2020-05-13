package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestUtil.testAnalysedReportData;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModelTestFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ReportableVariantAnalyzerTest {

    private static final double EPSILON = 1.0e-10;
    private static final String ONCO = "GENE1";
    private static final String TSG = "GENE2";

    private static final DriverGeneView TEST_DRIVER_GENE_VIEW = PatientReporterTestFactory.createTestDriverGeneView(ONCO, TSG);
    private static final List<DriverCatalog> TEST_DRIVER_CATALOG =
            Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalogEntry(ONCO),
                    PatientReporterTestFactory.createTestDriverCatalogEntry(TSG));

    @Test
    public void overruleDriverLikelihoodOfSomaticVariants() {
        SomaticVariant variant1 = PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(ONCO).build();
        SomaticVariant variant2 = PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build();
        List<SomaticVariant> variantsToReport = Lists.newArrayList(variant1, variant2);

        List<ReportableGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();

        ReportVariantAnalysis reportableVariantsAnalysis = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                TEST_DRIVER_CATALOG,
                TEST_DRIVER_GENE_VIEW,
                germlineVariantsToReport,
                PatientReporterTestFactory.createTestEmptyGermlineGenesReporting(),
                LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                testAnalysedReportData().actionabilityAnalyzer(),
                null);

        List<ReportableVariant> reportableVariants = reportableVariantsAnalysis.variantsToReport();

        assertEquals(4, reportableVariants.size());
        assertNullableEquals(1.0, reportableVariants.get(0).driverLikelihood(), EPSILON);
        assertNullableEquals(0.5, reportableVariants.get(1).driverLikelihood(), EPSILON);
        assertNullableEquals(1.0, reportableVariants.get(2).driverLikelihood(), EPSILON);
        assertNullableEquals(0.5, reportableVariants.get(3).driverLikelihood(), EPSILON);
    }

    @Test
    public void mergeWithoutGermlineVariantsWithDRUPActionabilityWorks() {
        SomaticVariant variant1 = PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(ONCO).build();
        SomaticVariant variant2 = PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build();

        List<SomaticVariant> variantsToReport = Lists.newArrayList(variant1, variant2);

        ReportVariantAnalysis reportableVariantsAnalysis = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                TEST_DRIVER_CATALOG,
                TEST_DRIVER_GENE_VIEW,
                Lists.newArrayList(),
                PatientReporterTestFactory.createTestEmptyGermlineGenesReporting(),
                LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                testAnalysedReportData().actionabilityAnalyzer(),
                null);

        List<ReportableVariant> reportableVariants = reportableVariantsAnalysis.variantsToReport();

        assertEquals(2, reportableVariants.size());

        assertEquals(ONCO, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());

        assertEquals(TSG, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithNotifyGermline() {
        List<SomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build());

        List<ReportableGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        ReportVariantAnalysis reportableVariantsAnalysis = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                TEST_DRIVER_CATALOG,
                TEST_DRIVER_GENE_VIEW,
                germlineVariantsToReport,
                germlineReportingModel,
                LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                testAnalysedReportData().actionabilityAnalyzer(),
                null);

        List<ReportableVariant> reportableVariants = reportableVariantsAnalysis.variantsToReport();

        assertEquals(3, reportableVariants.size());
        assertEquals(TSG, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());

        assertEquals(ONCO, reportableVariants.get(1).gene());
        assertTrue(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(TSG, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertFalse(reportableVariants.get(2).biallelic());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithoutNotifyGermline() {
        List<SomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build());

        List<ReportableGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        ReportVariantAnalysis reportableVariantsAnalysis = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                TEST_DRIVER_CATALOG,
                TEST_DRIVER_GENE_VIEW,
                germlineVariantsToReport,
                germlineReportingModel,
                LimsGermlineReportingLevel.NO_REPORTING,
                testAnalysedReportData().actionabilityAnalyzer(),
                null);

        List<ReportableVariant> reportableVariants = reportableVariantsAnalysis.variantsToReport();

        assertEquals(3, reportableVariants.size());
        assertEquals(TSG, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).biallelic());

        assertEquals(ONCO, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(TSG, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertFalse(reportableVariants.get(2).biallelic());
    }

    @NotNull
    private static GermlineReportingModel createGermlineReportingModelWithOncoAndTSG() {
        Map<String, Boolean> germlineGenesReportingMap = Maps.newHashMap();
        germlineGenesReportingMap.put(ONCO, true);
        germlineGenesReportingMap.put(TSG, false);
        return GermlineReportingModelTestFactory.buildFromMap(germlineGenesReportingMap);
    }

    @NotNull
    private static List<ReportableGermlineVariant> createBiallelicGermlineVariantsOnOncoAndTSG() {
        ReportableGermlineVariant germlineVariant1 = ImmutableReportableGermlineVariant.builder()
                .variant(PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(ONCO).biallelic(true).build())
                .driverLikelihood(1.0)
                .build();
        ReportableGermlineVariant germlineVariant2 = ImmutableReportableGermlineVariant.builder()
                .variant(PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(TSG).biallelic(false).build())
                .driverLikelihood(0.5)
                .build();
        return Lists.newArrayList(germlineVariant1, germlineVariant2);
    }

    private void assertNullableEquals(double expected, @Nullable Double actual, double epsilon) {
        assertNotNull(actual);
        assertEquals(expected, actual, epsilon);
    }
}