package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientreporter.algo.AnalysedReportData;
import com.hartwig.hmftools.patientreporter.algo.ImmutableAnalysedReportData;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.TaxonomyDb;
import com.hartwig.hmftools.patientreporter.virusbreakend.TaxonomyDbFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusInterpretationFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusInterpretationModel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestFactory {

    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PIPELINE_VERSION_FILE = RUN_DIRECTORY + "/pipeline.version";
    private static final String PURPLE_PURITY_TSV = RUN_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = RUN_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/sample.driver.catalog.somatic.tsv";
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/purple/sample.driver.catalog.germline.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = RUN_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String PURPLE_GERMLINE_VARIANT_VCF = RUN_DIRECTORY + "/purple/sample.purple.germline.vcf";
    private static final String PURPLE_SOMATIC_COPYNUMBER_TSV = RUN_DIRECTORY + "/purple/sample.purple.cnv.somatic.tsv";
    private static final String PURPLE_CIRCOS_FILE = RUN_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String LINX_FUSIONS_TSV = RUN_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = RUN_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/linx/sample.linx.driver.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = RUN_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = RUN_DIRECTORY + "/cuppa/sample.cuppa.conclusion.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_PLOT = RUN_DIRECTORY + "/cuppa/sample.cuppa.chart.png";
    private static final String VIRUS_BREAKEND_TSV = RUN_DIRECTORY + "/virusbreakend/sample.virusbreakend.vcf.summary.tsv";
    private static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/sample.peach.genotype.tsv";
    private static final String PROTECT_EVIDENCE_TSV = RUN_DIRECTORY + "/protect/sample.protect.tsv";

    private static final String SIGNATURE_PATH = Resources.getResource("signature/signature_test.png").getPath();
    private static final String RVA_LOGO_PATH = Resources.getResource("rva_logo/rva_logo_test.jpg").getPath();
    private static final String COMPANY_LOGO_PATH = Resources.getResource("company_logo/hartwig_logo_test.jpg").getPath();

    private static final String SAMPLE_SUMMARY_TSV = Resources.getResource("sample_summary/sample_summary.tsv").getPath();
    private static final String GERMLINE_REPORTING_TSV = Resources.getResource("germline_reporting/germline_reporting.tsv").getPath();
    private static final String TAXONOMY_DB_TSV = Resources.getResource("viral_reporting/taxonomy_db.tsv").getPath();
    private static final String VIRUS_INTERPRETATION_TSV = Resources.getResource("viral_reporting/virus_interpretation.tsv").getPath();
    private static final String VIRUS_BLACKLIST_TSV = Resources.getResource("viral_reporting/virus_blacklist.tsv").getPath();

    private PatientReporterTestFactory() {
    }

    @NotNull
    public static PatientReporterConfig createTestReporterConfig() {
        return ImmutablePatientReporterConfig.builder()
                .tumorSampleId(Strings.EMPTY)
                .tumorSampleBarcode(Strings.EMPTY)
                .outputDirReport(Strings.EMPTY)
                .outputDirData(Strings.EMPTY)
                .primaryTumorTsv(Strings.EMPTY)
                .limsDir(Strings.EMPTY)
                .rvaLogo(RVA_LOGO_PATH)
                .companyLogo(COMPANY_LOGO_PATH)
                .signature(SIGNATURE_PATH)
                .qcFail(false)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .purplePurityTsv(PURPLE_PURITY_TSV)
                .purpleQcFile(PURPLE_QC_FILE)
                .purpleSomaticDriverCatalogTsv(PURPLE_SOMATIC_DRIVER_CATALOG_TSV)
                .purpleGermlineDriverCatalogTsv(PURPLE_GERMLINE_DRIVER_CATALOG_TSV)
                .purpleSomaticVariantVcf(PURPLE_SOMATIC_VARIANT_VCF)
                .purpleGermlineVariantVcf(PURPLE_GERMLINE_VARIANT_VCF)
                .purpleSomaticCopyNumberTsv(PURPLE_SOMATIC_COPYNUMBER_TSV)
                .purpleCircosPlot(PURPLE_CIRCOS_FILE)
                .linxFusionTsv(LINX_FUSIONS_TSV)
                .linxBreakendTsv(LINX_BREAKEND_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .molecularTissueOriginTxt(MOLECULAR_TISSUE_ORIGIN_TXT)
                .molecularTissueOriginPlot(MOLECULAR_TISSUE_ORIGIN_PLOT)
                .virusBreakendTsv(VIRUS_BREAKEND_TSV)
                .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                .protectEvidenceTsv(PROTECT_EVIDENCE_TSV)
                .germlineReportingTsv(Strings.EMPTY)
                .sampleSummaryTsv(SAMPLE_SUMMARY_TSV)
                .taxonomyDbTsv(TAXONOMY_DB_TSV)
                .virusInterpretationTsv(VIRUS_INTERPRETATION_TSV)
                .virusBlacklistTsv(VIRUS_BLACKLIST_TSV)
                .isCorrectedReport(false)
                .onlyCreatePDF(false)
                .expectedPipelineVersion("5.22")
                .overridePipelineVersion(false)
                .build();
    }

    @NotNull
    public static ReportData loadTestReportData() {
        List<PatientPrimaryTumor> patientPrimaryTumors = Lists.newArrayList();
        Lims lims = LimsFactory.empty();

        return ImmutableQCFailReportData.builder()
                .patientPrimaryTumors(patientPrimaryTumors)
                .limsModel(lims)
                .signaturePath(SIGNATURE_PATH)
                .logoRVAPath(RVA_LOGO_PATH)
                .logoCompanyPath(COMPANY_LOGO_PATH)
                .build();
    }

    @NotNull
    public static AnalysedReportData loadTestAnalysedReportData() {
        try {
            GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(GERMLINE_REPORTING_TSV);
            SummaryModel summaryModel = SummaryFile.buildFromTsv(SAMPLE_SUMMARY_TSV);
            TaxonomyDb taxonomyDb = TaxonomyDbFile.loadFromTsv(TAXONOMY_DB_TSV);
            VirusInterpretationModel virusInterpretationModel = VirusInterpretationFile.buildFromTsv(VIRUS_INTERPRETATION_TSV);
            VirusBlacklistModel virusBlackListModel = VirusBlacklistFile.buildFromTsv(VIRUS_BLACKLIST_TSV);

            return ImmutableAnalysedReportData.builder()
                    .from(loadTestReportData())
                    .germlineReportingModel(germlineReportingModel)
                    .summaryModel(summaryModel)
                    .taxonomyDb(taxonomyDb)
                    .virusInterpretationModel(virusInterpretationModel)
                    .virusBlackListModel(virusBlackListModel)
                    .build();
        } catch (IOException exception) {
            throw new IllegalStateException("Could not load test analysed report data: " + exception.getMessage());
        }
    }
}
