package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModelTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.germline.ImmutableInterpretGermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.germline.InterpretGermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantAnalyzerTest {

    private static final String GENE_1 = "GENE1";
    private static final String GENE_2 = "GENE2";

    private static final DriverGeneView TEST_DRIVER_GENE_VIEW =
            PatientReporterTestFactory.createTestDriverGeneView(GENE_1, GENE_2);
    private static final List<DriverCatalog> TEST_DRIVER_CATALOG =
            Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalogEntry(GENE_1),
                    PatientReporterTestFactory.createTestDriverCatalogEntry(GENE_2));

    @Test
    public void mergeWithoutGermlineVariantsWithDRUPActionabilityWorks() {
        SomaticVariant variant1 = PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_1).build();
        SomaticVariant variant2 = PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_2).build();

        List<SomaticVariant> variantsToReport = Lists.newArrayList(variant1, variant2);

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                TEST_DRIVER_CATALOG,
                TEST_DRIVER_GENE_VIEW,
                Lists.newArrayList(),
                PatientReporterTestFactory.createTestEmptyGermlineGenesReporting(),
                LimsGermlineReportingChoice.ALL);

        assertEquals(2, reportableVariants.size());

        assertEquals(GENE_1, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());

        assertEquals(GENE_2, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithNotifyGermline() {
        List<SomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_2).build());

        List<InterpretGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                TEST_DRIVER_CATALOG,
                TEST_DRIVER_GENE_VIEW,
                germlineVariantsToReport,
                germlineReportingModel,
                LimsGermlineReportingChoice.ALL);

        assertEquals(3, reportableVariants.size());
        assertEquals(GENE_2, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());

        assertEquals(GENE_1, reportableVariants.get(1).gene());
        assertTrue(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(GENE_2, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(2).biallelic());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithoutNotifyGermline() {
        List<SomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_2).build());

        List<InterpretGermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                TEST_DRIVER_CATALOG,
                TEST_DRIVER_GENE_VIEW,
                germlineVariantsToReport,
                germlineReportingModel,
                LimsGermlineReportingChoice.NONE_ALLOW_FAMILY);

        assertEquals(3, reportableVariants.size());
        assertEquals(GENE_2, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).biallelic());

        assertEquals(GENE_1, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(GENE_2, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(2).biallelic());
    }

    @NotNull
    private static GermlineReportingModel createGermlineReportingModelWithOncoAndTSG() {
        Map<String, Boolean> germlineGenesReportingMap = Maps.newHashMap();
        germlineGenesReportingMap.put(GENE_1, true);
        germlineGenesReportingMap.put(GENE_2, false);
        return GermlineReportingModelTestFactory.buildFromMap(germlineGenesReportingMap);
    }

    @NotNull
    private static List<InterpretGermlineVariant> createBiallelicGermlineVariantsOnOncoAndTSG() {
        InterpretGermlineVariant germlineVariant1 = ImmutableInterpretGermlineVariant.builder().germlineVariant(
                PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(GENE_1).biallelic(true).build()).driverLikelihood(0).build();
        InterpretGermlineVariant germlineVariant2 =ImmutableInterpretGermlineVariant.builder().germlineVariant(
                PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(GENE_2).biallelic(true).build()).driverLikelihood(0).build();
        return Lists.newArrayList(germlineVariant1, germlineVariant2);
    }
}