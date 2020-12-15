package com.hartwig.hmftools.patientreporter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.WithEvent;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.protect.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class AnalysedPatientReporterTest {

    private static final String TUMOR_SAMPLE_ID = "sample";
    private static final String REF_SAMPLE_ID = "ref_sample";
    private static final String PATIENT_ID = "patient";

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = BASE_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_DRIVER_CATALOG_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String BACHELOR_TSV = BASE_DIRECTORY + "/bachelor/sample.reportable_germline_variants.tsv";
    private static final String LINX_FUSIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = BASE_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_VIRAL_INSERTION_TSV = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";
    private static final String LINX_DRIVERS_TSV = BASE_DIRECTORY + "/linx/sample.drivers.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String CIRCOS_FILE = BASE_DIRECTORY + "/purple/plot/sample.circos.png";

    @Test
    public void canRunOnRunDirectory() throws IOException {
        AnalysedPatientReporter reporter = new AnalysedPatientReporter(PatientReporterTestFactory.loadTestAnalysedReportData());

        SampleMetadata sampleMetadata = ImmutableSampleMetadata.builder()
                .patientId(PATIENT_ID)
                .refSampleId(REF_SAMPLE_ID)
                .refSampleBarcode(Strings.EMPTY)
                .tumorSampleId(TUMOR_SAMPLE_ID)
                .tumorSampleBarcode(Strings.EMPTY)
                .build();

        assertNotNull(reporter.run(sampleMetadata,
                PURPLE_PURITY_TSV,
                PURPLE_QC_FILE,
                PURPLE_DRIVER_CATALOG_TSV,
                PURPLE_SOMATIC_VARIANT_VCF,
                BACHELOR_TSV,
                LINX_FUSIONS_TSV,
                LINX_BREAKEND_TSV,
                LINX_VIRAL_INSERTION_TSV,
                LINX_DRIVERS_TSV,
                CHORD_PREDICTION_TXT,
                CIRCOS_FILE,
                null,
                false));
    }

    @Test
    public void canDetermineForNumber() {
        double purityCorrect = 0.40;
        boolean hasReliablePurityCorrect = true;

        assertEquals(QsFormNumber.FOR_080.display(), AnalysedPatientReporter.determineForNumber(hasReliablePurityCorrect, purityCorrect));

        double purityNotCorrect = 0.10;
        boolean hasReliablePurityNotCorrect = false;

        assertEquals(QsFormNumber.FOR_209.display(),
                AnalysedPatientReporter.determineForNumber(hasReliablePurityNotCorrect, purityNotCorrect));
    }

    @Test
    public void canFilterVariantsForGermlineConsent() {
        ReportableVariant somaticVariant = createTestReportableVariantBuilder().source(ReportableVariantSource.SOMATIC).build();
        ReportableVariant germlineVariant = createTestReportableVariantBuilder().source(ReportableVariantSource.GERMLINE).build();

        assertEquals(2,
                AnalysedPatientReporter.filterVariantsForGermlineConsent(Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION).size());

        assertEquals(1,
                AnalysedPatientReporter.filterVariantsForGermlineConsent(Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.NO_REPORTING).size());
    }

    @Test
    public void canFilterEvidenceForGermlineConsent() {
        String somaticEventString = "somatic";
        String germlineEventString = "germline";

        ReportableVariant somaticVariant = variant(ReportableVariantSource.SOMATIC, somaticEventString);
        ReportableVariant germlineVariant = variant(ReportableVariantSource.GERMLINE, germlineEventString);

        // Extra space is added because of protein impact which is missing in this test.
        WithEvent somaticEvent = new AnalysedPatientReporterTest.DummyEvent(somaticEventString + " ");
        WithEvent germlineEvent = new AnalysedPatientReporterTest.DummyEvent(germlineEventString + " ");

        assertEquals(2,
                AnalysedPatientReporter.filterEvidenceForGermlineConsent(Lists.newArrayList(somaticEvent, germlineEvent),
                        Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION).size());

        assertEquals(1,
                AnalysedPatientReporter.filterEvidenceForGermlineConsent(Lists.newArrayList(somaticEvent, germlineEvent),
                        Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.NO_REPORTING).size());

        assertEquals(2,
                AnalysedPatientReporter.filterEvidenceForGermlineConsent(Lists.newArrayList(somaticEvent, somaticEvent),
                        Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.NO_REPORTING).size());

        assertEquals(2,
                AnalysedPatientReporter.filterEvidenceForGermlineConsent(Lists.newArrayList(somaticEvent, germlineEvent),
                        Lists.newArrayList(somaticVariant),
                        LimsGermlineReportingLevel.NO_REPORTING).size());
    }

    @NotNull
    private static ReportableVariant variant(@NotNull ReportableVariantSource source, @NotNull String event) {
        return createTestReportableVariantBuilder().source(source).gene(event).build();
    }

    @NotNull
    private static ImmutableReportableVariant.Builder createTestReportableVariantBuilder() {
        return ImmutableReportableVariant.builder()
                .source(ReportableVariantSource.SOMATIC)
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .totalCopyNumber(0)
                .alleleCopyNumber(0D)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0D)
                .biallelic(false);
    }

    private static class DummyEvent implements WithEvent {

        @NotNull
        private final String event;

        public DummyEvent(@NotNull final String event) {
            this.event = event;
        }

        @NotNull
        @Override
        public String event() {
            return event;
        }
    }
}
