package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantAnalyzerTest {

    private static final String ONCOGENE = "ONCO";
    private static final String TSG = "TSG";

    @Test
    public void mergeWithoutGermlineVariant() {
        List<EnrichedSomaticVariant> variantsToReport = Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().build());
        List<DriverCatalog> driverCatalog = Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalog().build());
        Map<String, DriverCategory> driverCategoryMap = PatientReporterTestFactory.createTestDriverCategoryMap();
        Set<String> drupActionableGenes = Sets.newHashSet(ONCOGENE, TSG);

        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variantsToReport) {

            DriverCatalog catalog = ReportableVariantAnalyzer.catalogEntryForVariant(driverCatalog, variant.gene());

            reportableVariants.add(ReportableVariantAnalyzer.fromSomaticVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryMap.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .notifyClinicalGeneticist(false)
                    .build());
        }
        assertEquals(1, reportableVariants.size());

    }

    @Test
    public void mergeSomaticAndGermlineVariant() {
        List<EnrichedSomaticVariant> variantsToReport = Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().build());
        List<DriverCatalog> driverCatalog = Lists.newArrayList(PatientReporterTestFactory.createTestDriverCatalog().build());
        Map<String, DriverCategory> driverCategoryMap = PatientReporterTestFactory.createTestDriverCategoryMap();
        Set<String> drupActionableGenes = Sets.newHashSet(ONCOGENE, TSG);
        List<GermlineVariant > filteredGermlineVariants =  createTestGermlineVariantsONCOGene();
        Map<String, Boolean> germlineGenesToNotifyMap = PatientReporterTestFactory.createTestGermlineGenesReporting();


        List<ReportableVariant> reportableVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : variantsToReport) {

            DriverCatalog catalog = ReportableVariantAnalyzer.catalogEntryForVariant(driverCatalog, variant.gene());


            reportableVariants.add(ReportableVariantAnalyzer.fromSomaticVariant(variant).isDrupActionable(drupActionableGenes.contains(variant.gene()))
                    .driverCategory(driverCategoryMap.get(variant.gene()))
                    .driverLikelihood(catalog != null ? catalog.driverLikelihood() : null)
                    .notifyClinicalGeneticist(false)
                    .build());
        }

        for (GermlineVariant germlineVariant : filteredGermlineVariants) {
            reportableVariants.add(ReportableVariantAnalyzer.fromGermlineVariant(germlineVariant).isDrupActionable(drupActionableGenes.contains(germlineVariant.gene()))
                    .driverCategory(driverCategoryMap.get(germlineVariant.gene()))
                    .driverLikelihood(null)
                    .notifyClinicalGeneticist(germlineGenesToNotifyMap.get(germlineVariant.gene()))
                    .build());

        }

        assertEquals(2, reportableVariants.size());

    }


    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsONCOGene() {
        return Lists.newArrayList(PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(ONCOGENE).build());
    }

}