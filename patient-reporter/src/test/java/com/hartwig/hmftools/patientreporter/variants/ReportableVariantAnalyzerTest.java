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
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.DriverGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModelTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.ImmutableDriverGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.somatic.DriverSomaticVariant;
import com.hartwig.hmftools.patientreporter.variants.somatic.ImmutableDriverSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ReportableVariantAnalyzerTest {

    private static final double EPSILON = 1.0e-10;
    private static final String ONCO = "KIT";
    private static final String TSG = "PTEN";

    private static final List<DriverCatalog> TEST_DRIVER_CATALOG =
            Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalogEntry(ONCO),
                    PatientReporterTestFactory.createTestDriverCatalogEntry(TSG));

    @Test
    public void overruleDriverLikelihoodOfSomaticVariants() {
        DriverSomaticVariant variant1 = ImmutableDriverSomaticVariant.builder()
                .variant(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(ONCO).build())
                .driverLikelihood(0D)
                .build();
        DriverSomaticVariant variant2 = ImmutableDriverSomaticVariant.builder()
                .variant(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build())
                .driverLikelihood(0D)
                .build();
        List<DriverSomaticVariant> somaticVariantsToReport = Lists.newArrayList(variant1, variant2);
        List<DriverGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();

        ReportableVariantAnalysis reportableVariantsAnalysis =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantsToReport,
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
        DriverSomaticVariant variant1 = ImmutableDriverSomaticVariant.builder()
                .variant(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(ONCO).build())
                .driverLikelihood(0D)
                .build();
        DriverSomaticVariant variant2 = ImmutableDriverSomaticVariant.builder()
                .variant(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build())
                .driverLikelihood(0D)
                .build();

        List<DriverSomaticVariant> somaticVariantsToReport = Lists.newArrayList(variant1, variant2);

        ReportableVariantAnalysis reportableVariantsAnalysis =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantsToReport,
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
        List<DriverSomaticVariant> somaticVariantsToReport = Lists.newArrayList(ImmutableDriverSomaticVariant.builder()
                .variant(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build())
                .driverLikelihood(0D)
                .build());

        List<DriverGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        ReportableVariantAnalysis reportableVariantsAnalysis =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantsToReport,
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
        List<DriverSomaticVariant> somaticVariantsToReport = Lists.newArrayList(ImmutableDriverSomaticVariant.builder()
                .variant(PatientReporterTestFactory.createTestSomaticVariantBuilder().gene(TSG).build())
                .driverLikelihood(0D)
                .build());

        List<DriverGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        ReportableVariantAnalysis reportableVariantsAnalysis =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantsToReport,
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
    private static List<DriverGermlineVariant> createBiallelicGermlineVariantsOnOncoAndTSG() {
        DriverGermlineVariant germlineVariant1 = ImmutableDriverGermlineVariant.builder()
                .variant(PatientReporterTestFactory.createTestGermlineVariantBuilder()
                        .gene(ONCO)
                        .biallelic(true)
                        .adjustedVaf(0)
                        .adjustedCopyNumber(0)
                        .build())
                .driverLikelihood(1.0)
                .build();
        DriverGermlineVariant germlineVariant2 = ImmutableDriverGermlineVariant.builder()
                .variant(PatientReporterTestFactory.createTestGermlineVariantBuilder()
                        .gene(TSG)
                        .biallelic(false)
                        .adjustedVaf(0)
                        .adjustedCopyNumber(0)
                        .build())
                .driverLikelihood(0.5)
                .build();
        return Lists.newArrayList(germlineVariant1, germlineVariant2);
    }

    private void assertNullableEquals(double expected, @Nullable Double actual, double epsilon) {
        assertNotNull(actual);
        assertEquals(expected, actual, epsilon);
    }
}