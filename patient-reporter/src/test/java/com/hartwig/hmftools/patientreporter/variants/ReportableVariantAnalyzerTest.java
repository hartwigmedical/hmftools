package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantAnalyzerTest {

    private static final String ONCOGENE = "ONCO";
    private static final String TSG = "TSG";

    @Test
    public void mergeWithoutGermlineVariant() {
        EnrichedSomaticVariant variants1 = PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(ONCOGENE).build();
        EnrichedSomaticVariant variants2 = PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(TSG).build();
        List<EnrichedSomaticVariant> variantsToReport = Lists.newArrayList(variants1, variants2);
        List<DriverCatalog> driverCatalog = Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalog().build());
        Map<String, DriverCategory> driverCategoryMap = PatientReporterTestFactory.createTestDriverCategoryMap();
        Set<String> drupActionableGenes = Sets.newHashSet(ONCOGENE, TSG);

        GermlineReportingModel germlineGenesreportingModel = PatientReporterTestFactory.createTestEmptyGermlineGenesReporting();


        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                driverCatalog,
                driverCategoryMap,
                drupActionableGenes,
                Lists.newArrayList(),
                germlineGenesreportingModel, LimsGermlineReportingChoice.ALL);

        assertEquals(2, reportableVariants.size());
        assertEquals(ONCOGENE, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).biallelic());

        assertEquals(TSG, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithNotifyGermline() {
        List<EnrichedSomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(TSG).build());
        List<DriverCatalog> driverCatalog = Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalog().build());
        Map<String, DriverCategory> driverCategoryMap = PatientReporterTestFactory.createTestDriverCategoryMap();
        Set<String> drupActionableGenes = Sets.newHashSet(ONCOGENE, TSG);

        List<GermlineVariant> filteredGermlineVariants = createTestGermlineVariantsONCOGene();
        GermlineReportingModel germlineGenesreportingModel = PatientReporterTestFactory.createTestGermlineGenesReporting();

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                driverCatalog,
                driverCategoryMap,
                drupActionableGenes,
                filteredGermlineVariants,
                germlineGenesreportingModel, LimsGermlineReportingChoice.ALL);

        assertEquals(3, reportableVariants.size());
        assertEquals(TSG, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).biallelic());

        assertEquals(ONCOGENE, reportableVariants.get(1).gene());
        assertTrue(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(TSG, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(2).biallelic());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithoutNotifyGermline() {
        List<EnrichedSomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(TSG).build());
        List<DriverCatalog> driverCatalog = Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalog().build());
        Map<String, DriverCategory> driverCategoryMap = PatientReporterTestFactory.createTestDriverCategoryMap();
        Set<String> drupActionableGenes = Sets.newHashSet(ONCOGENE, TSG);

        List<GermlineVariant> filteredGermlineVariants = createTestGermlineVariantsONCOGene();
        GermlineReportingModel germlineGenesreportingModel = PatientReporterTestFactory.createTestGermlineGenesReporting();

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                driverCatalog,
                driverCategoryMap,
                drupActionableGenes,
                filteredGermlineVariants,
                germlineGenesreportingModel, LimsGermlineReportingChoice.NONE);

        assertEquals(3, reportableVariants.size());
        assertEquals(TSG, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).biallelic());

        assertEquals(ONCOGENE, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(TSG, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(2).biallelic());
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsONCOGene() {
        GermlineVariant germlineVariant1 =
                PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(ONCOGENE).biallelic(true).build();
        GermlineVariant germlineVariant2 = PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(TSG).biallelic(true).build();
        return Lists.newArrayList(germlineVariant1, germlineVariant2);
    }

}