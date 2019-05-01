package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModelTestFactory;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableVariantAnalyzerTest {

    private static final String GENE_WITH_NOTIFY = "GENE1";
    private static final String GENE_WITHOUT_NOTIFY = "GENE2";

    @Test
    public void mergeWithoutGermlineVariantsWithDRUPActionabilityWorks() {
        EnrichedSomaticVariant variant1 =
                PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_WITH_NOTIFY).build();
        EnrichedSomaticVariant variant2 =
                PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_WITHOUT_NOTIFY).build();

        List<EnrichedSomaticVariant> variantsToReport = Lists.newArrayList(variant1, variant2);
        Set<String> drupActionableGenes = Sets.newHashSet(GENE_WITH_NOTIFY);

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                Lists.newArrayList(),
                Maps.newHashMap(),
                drupActionableGenes,
                Lists.newArrayList(),
                PatientReporterTestFactory.createTestEmptyGermlineGenesReporting(),
                LimsGermlineReportingChoice.ALL);

        assertEquals(2, reportableVariants.size());

        assertEquals(GENE_WITH_NOTIFY, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).isDrupActionable());

        assertEquals(GENE_WITHOUT_NOTIFY, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
        assertFalse(reportableVariants.get(1).isDrupActionable());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithNotifyGermline() {
        List<EnrichedSomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_WITHOUT_NOTIFY).build());

        List<GermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                Lists.newArrayList(),
                Maps.newHashMap(),
                Sets.newHashSet(),
                germlineVariantsToReport,
                germlineReportingModel,
                LimsGermlineReportingChoice.ALL);

        assertEquals(3, reportableVariants.size());
        assertEquals(GENE_WITHOUT_NOTIFY, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());

        assertEquals(GENE_WITH_NOTIFY, reportableVariants.get(1).gene());
        assertTrue(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(GENE_WITHOUT_NOTIFY, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(2).biallelic());
    }

    @Test
    public void mergeSomaticAndGermlineVariantWithoutNotifyGermline() {
        List<EnrichedSomaticVariant> variantsToReport =
                Lists.newArrayList(PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene(GENE_WITHOUT_NOTIFY).build());

        List<GermlineVariant> germlineVariantsToReport = createBiallelicGermlineVariantsOnOncoAndTSG();
        GermlineReportingModel germlineReportingModel = createGermlineReportingModelWithOncoAndTSG();

        List<ReportableVariant> reportableVariants = ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(variantsToReport,
                Lists.newArrayList(),
                Maps.newHashMap(),
                Sets.newHashSet(),
                germlineVariantsToReport,
                germlineReportingModel,
                LimsGermlineReportingChoice.NONE_ALLOW_FAMILY);

        assertEquals(3, reportableVariants.size());
        assertEquals(GENE_WITHOUT_NOTIFY, reportableVariants.get(0).gene());
        assertFalse(reportableVariants.get(0).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(0).biallelic());

        assertEquals(GENE_WITH_NOTIFY, reportableVariants.get(1).gene());
        assertFalse(reportableVariants.get(1).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(1).biallelic());

        assertEquals(GENE_WITHOUT_NOTIFY, reportableVariants.get(2).gene());
        assertFalse(reportableVariants.get(2).notifyClinicalGeneticist());
        assertTrue(reportableVariants.get(2).biallelic());
    }

    @NotNull
    private static GermlineReportingModel createGermlineReportingModelWithOncoAndTSG() {
        Map<String, Boolean> germlineGenesReportingMap = Maps.newHashMap();
        germlineGenesReportingMap.put(GENE_WITH_NOTIFY, true);
        germlineGenesReportingMap.put(GENE_WITHOUT_NOTIFY, false);
        return GermlineReportingModelTestFactory.buildFromMap(germlineGenesReportingMap);
    }

    @NotNull
    private static List<GermlineVariant> createBiallelicGermlineVariantsOnOncoAndTSG() {
        GermlineVariant germlineVariant1 =
                PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(GENE_WITH_NOTIFY).biallelic(true).build();
        GermlineVariant germlineVariant2 =
                PatientReporterTestFactory.createTestGermlineVariantBuilder().gene(GENE_WITHOUT_NOTIFY).biallelic(true).build();
        return Lists.newArrayList(germlineVariant1, germlineVariant2);
    }
}