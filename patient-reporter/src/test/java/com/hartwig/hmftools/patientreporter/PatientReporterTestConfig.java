package com.hartwig.hmftools.patientreporter;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestConfig {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PIPELINE_VERSION_FILE = BASE_DIRECTORY + "/pipeline.version";
    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = BASE_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_SOMATIC_DRIVER_CATALOG_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.somatic.tsv";
    private static final String PURPLE_GERMLINE_DRIVER_CATALOG_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.germline.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String PURPLE_GERMLINE_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.germline.vcf";
    private static final String PURPLE_SOMATIC_COPYNUMBER_TSV = BASE_DIRECTORY + "/purple/sample.purple.cnv.somatic.tsv";
    private static final String PURPLE_CIRCOS_FILE = BASE_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String LINX_FUSIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = BASE_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_DRIVER_CATALOG_TSV = BASE_DIRECTORY + "/linx/sample.linx.driver.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = BASE_DIRECTORY + "/cuppa/sample.cuppa.conclusion.txt";
    private static final String MOLECULAR_TISSUE_ORIGIN_PLOT = BASE_DIRECTORY + "/purple/plot/sample.circos.png";
    private static final String VIRUS_BREAKEND_TSV = BASE_DIRECTORY + "/virusbreakend/sample.virusbreakend.vcf.summary.tsv";
    private static final String PEACH_GENOTYPE_TSV = BASE_DIRECTORY + "/peach/sample.peach.genotype.tsv";
    private static final String PROTECT_EVIDENCE_TSV = BASE_DIRECTORY + "/protect/sample.protect.tsv";

    private static final String VIRUS_DB_TSV = Resources.getResource("virusbreakend/virusdb.tsv").getPath();
    private static final String VIRUS_SUMMARY_TSV = Resources.getResource("virusbreakend/virus_summary.tsv").getPath();
    private static final String VIRUS_BLACKLIST_TSV = Resources.getResource("virusbreakend/virus_blacklist.tsv").getPath();

    private PatientReporterTestConfig() {
    }

    @NotNull
    public static PatientReporterConfig create() {
        return ImmutablePatientReporterConfig.builder()
                .tumorSampleId(Strings.EMPTY)
                .tumorSampleBarcode(Strings.EMPTY)
                .outputDirReport(Strings.EMPTY)
                .outputDirData(Strings.EMPTY)
                .reportingDbTsv(Strings.EMPTY)
                .primaryTumorTsv(Strings.EMPTY)
                .limsDir(Strings.EMPTY)
                .rvaLogo(Strings.EMPTY)
                .companyLogo(Strings.EMPTY)
                .signature(Strings.EMPTY)
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
                .sampleSummaryTsv(Strings.EMPTY)
                .virusDbTsv(VIRUS_DB_TSV)
                .virusSummaryTsv(VIRUS_SUMMARY_TSV)
                .virusBlacklistTsv(VIRUS_BLACKLIST_TSV)
                .isCorrectedReport(false)
                .onlyCreatePDF(false)
                .build();
    }
}
