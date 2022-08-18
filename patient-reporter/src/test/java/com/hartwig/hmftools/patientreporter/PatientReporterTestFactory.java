package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientreporter.algo.AnalysedReportData;
import com.hartwig.hmftools.patientreporter.algo.ImmutableAnalysedReportData;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReportData;
import com.hartwig.hmftools.patientreporter.remarks.SpecialRemarkFile;
import com.hartwig.hmftools.patientreporter.remarks.SpecialRemarkModel;

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
    private static final String LINX_SVS_TSV = RUN_DIRECTORY + "/linx/sample.linx.svs.tsv";
    private static final String LINX_DRIVER_CATALOG_TSV = RUN_DIRECTORY + "/linx/sample.linx.driver.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = RUN_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = RUN_DIRECTORY + "/cuppa/sample.cuppa.conclusion.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_PLOT = RUN_DIRECTORY + "/cuppa/sample.cuppa.chart.png";
    private static final String ANNOTATED_VIRUS_TSV = RUN_DIRECTORY + "/virusbreakend/sample.virus.annotated.tsv";
    private static final String PEACH_GENOTYPE_TSV = RUN_DIRECTORY + "/peach/sample.peach.genotype.tsv";
    private static final String PROTECT_EVIDENCE_TSV = RUN_DIRECTORY + "/protect/sample.protect.tsv";
    private static final String LILAC_RESULT_CSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.csv";
    private static final String LILAC_QC_CSV = RUN_DIRECTORY + "/lilac/tumor_sample.lilac.qc.csv";
    private static final String ROSE_TSV = RUN_DIRECTORY + "/rose/tumor_sample.rose.tsv";
    private static final String SIGNATURE_PATH = Resources.getResource("signature/signature_test.png").getPath();
    private static final String RVA_LOGO_PATH = Resources.getResource("rva_logo/rva_logo_test.jpg").getPath();
    private static final String COMPANY_LOGO_PATH = Resources.getResource("company_logo/hartwig_logo_test.jpg").getPath();

    private static final String SAMPLE_SPECIAL_REMARK_TSV = Resources.getResource("special_remark/sample_special_remark.tsv").getPath();

    private static final String GERMLINE_REPORTING_TSV = Resources.getResource("germline_reporting/germline_reporting.tsv").getPath();
    private static final String KNOWN_FUSION_FILE = Resources.getResource("known_fusion_data/known_fusion_file.csv").getPath();

    private static final String UDI_DI = "(01)8720299486027(8012)v5.28";

    private PatientReporterTestFactory() {
    }

    @NotNull
    public static PatientReporterConfig createTestReporterConfig() {
        return ImmutablePatientReporterConfig.builder()
                .tumorSampleId(Strings.EMPTY)
                .tumorSampleBarcode(Strings.EMPTY)
                .outputDirReport(Strings.EMPTY)
                .outputDirData(Strings.EMPTY)
                .knownFusionFile(KNOWN_FUSION_FILE)
                .primaryTumorTsv(Strings.EMPTY)
                .limsDir(Strings.EMPTY)
                .rvaLogo(RVA_LOGO_PATH)
                .companyLogo(COMPANY_LOGO_PATH)
                .signature(SIGNATURE_PATH)
                .udiDi(UDI_DI)
                .qcFail(false)
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
                .linxSvsTsv(LINX_SVS_TSV)
                .linxDriverCatalogTsv(LINX_DRIVER_CATALOG_TSV)
                .chordPredictionTxt(CHORD_PREDICTION_TXT)
                .molecularTissueOriginTxt(MOLECULAR_TISSUE_ORIGIN_TXT)
                .molecularTissueOriginPlot(MOLECULAR_TISSUE_ORIGIN_PLOT)
                .annotatedVirusTsv(ANNOTATED_VIRUS_TSV)
                .peachGenotypeTsv(PEACH_GENOTYPE_TSV)
                .protectEvidenceTsv(PROTECT_EVIDENCE_TSV)
                .lilacResultCsv(LILAC_RESULT_CSV)
                .lilacQcCsv(LILAC_QC_CSV)
                .addRose(true)
                .roseTsv(ROSE_TSV)
                .germlineReportingTsv(Strings.EMPTY)
                .sampleSpecialRemarkTsv(SAMPLE_SPECIAL_REMARK_TSV)
                .isCorrectedReport(false)
                .isCorrectedReportExtern(false)
                .onlyCreatePDF(false)
                .requirePipelineVersionFile(true)
                .pipelineVersionFile(PIPELINE_VERSION_FILE)
                .expectedPipelineVersion("5.28")
                .overridePipelineVersion(false)
                .refGenomeVersion(RefGenomeVersion.V37)
                .allowDefaultCohortConfig(false)
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
                .udiDi(UDI_DI)
                .build();
    }

    @NotNull
    public static AnalysedReportData loadTestAnalysedReportData() {
        try {
            GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(GERMLINE_REPORTING_TSV);
            SpecialRemarkModel specialRemarkModel = SpecialRemarkFile.buildFromTsv(SAMPLE_SPECIAL_REMARK_TSV);
            KnownFusionCache knownFusionCache = new KnownFusionCache();

            return ImmutableAnalysedReportData.builder()
                    .from(loadTestReportData())
                    .germlineReportingModel(germlineReportingModel)
                    .specialRemarkModel(specialRemarkModel)
                    .knownFusionCache(knownFusionCache)
                    .build();
        } catch (IOException exception) {
            throw new IllegalStateException("Could not load test analysed report data: " + exception.getMessage());
        }
    }
}
