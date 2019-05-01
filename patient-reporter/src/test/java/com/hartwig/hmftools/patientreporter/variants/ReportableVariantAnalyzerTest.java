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

        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variantsToReport) {
            DriverCatalog catalog = ReportableVariantAnalyzer.catalogEntryForVariant(driverCatalog, variant.gene());

            reportableVariants.add(ReportableVariantAnalyzer.fromSomaticVariant(variant)
                    .isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryMap.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        assertEquals(2, reportableVariants.size());
        assertEquals(ONCOGENE, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).biallelic());

        assertEquals(TSG, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());
    }

    @Test
    public void mergeSomaticAndGermlineVariant() {
        List<EnrichedSomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(TSG).build());
        List<DriverCatalog> driverCatalog = Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalog().build());
        Map<String, DriverCategory> driverCategoryMap = PatientReporterTestFactory.createTestDriverCategoryMap();
        Set<String> drupActionableGenes = Sets.newHashSet(ONCOGENE, TSG);

        List<GermlineVariant> filteredGermlineVariants = createTestGermlineVariantsONCOGene();
        GermlineReportingModel germlineGenesreportingModel = PatientReporterTestFactory.createTestGermlineGenesReporting();

        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variantsToReport) {

            DriverCatalog catalog = ReportableVariantAnalyzer.catalogEntryForVariant(driverCatalog, variant.gene());

            reportableVariants.add(ReportableVariantAnalyzer.fromSomaticVariant(variant)
                    .isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryMap.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        for (GermlineVariant germlineVariant : filteredGermlineVariants) {
            reportableVariants.add(ReportableVariantAnalyzer.fromGermlineVariant(germlineVariant)
                    .isDrupActionable(drupActionableGenes.contains(germlineVariant.gene()))
                    .driverCategory(driverCategoryMap.get(germlineVariant.gene()))
                    .driverLikelihood(null)
                    .notifyClinicalGeneticist(germlineGenesreportingModel.genesToNotifyClinicalGeneticist(germlineVariant.gene()))
                    .build());

        }

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

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsONCOGene() {
        GermlineVariant germlineVariant1 =
                PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(ONCOGENE).biallelic(true).build();
        GermlineVariant germlineVariant2 = PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(TSG).biallelic(true).build();
        return Lists.newArrayList(germlineVariant1, germlineVariant2);
    }

}