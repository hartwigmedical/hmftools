package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.nio.file.Files;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.patientreporter.panel.PanelFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PatientReporterConfig {

    Logger LOGGER = LogManager.getLogger(PatientReporterConfig.class);

    // General params needed for every report
    String TUMOR_SAMPLE_ID = "tumor_sample_id";
    String TUMOR_SAMPLE_BARCODE = "tumor_sample_barcode";
    String OUTPUT_DIRECTORY_REPORT = "output_dir_report";
    String OUTPUT_DIRECTORY_DATA = "output_dir_data";

    String PRIMARY_TUMOR_TSV = "primary_tumor_tsv";
    String LIMS_DIRECTORY = "lims_dir";
    String PEACH_GENOTYPE_TSV = "peach_genotype_tsv";

    String RVA_LOGO = "rva_logo";
    String COMPANY_LOGO = "company_logo";
    String SIGNATURE = "signature";

    String UDI_DI = "udi_di";

    // General params needed for every report but for QC fail it can be optional in some cases
    String REF_SAMPLE_ID = "ref_sample_id";
    String REF_SAMPLE_BARCODE = "ref_sample_barcode";

    // Params specific for QC Fail reports
    String QC_FAIL = "qc_fail";
    String QC_FAIL_REASON = "qc_fail_reason";

    // Params specific for Panel reports
    String PANEL = "panel";
    String PANEL_QC_FAIL = "panel_qc_fail";
    String PANEL_QC_FAIL_REASON = "panel_qc_fail_reason";
    String PANEL_VCF_NAME = "panel_vcf_name";

    // Params specific for actual patient reports
    String PURPLE_PURITY_TSV = "purple_purity_tsv"; // Also used for certain QC fail reports in case deep WGS is available.
    String PURPLE_QC_FILE = "purple_qc_file"; // Also used for certain QC fail reports in case deep WGS is available.
    String PURPLE_SOMATIC_DRIVER_CATALOG_TSV = "purple_somatic_driver_catalog_tsv";
    String PURPLE_GERMLINE_DRIVER_CATALOG_TSV = "purple_germline_driver_catalog_tsv";
    String PURPLE_SOMATIC_VARIANT_VCF = "purple_somatic_variant_vcf";
    String PURPLE_GERMLINE_VARIANT_VCF = "purple_germline_variant_vcf";
    String PURPLE_SOMATIC_COPYNUMBER_TSV = "purple_somatic_copynumber_tsv";
    String PURPLE_CIRCOS_PLOT = "purple_circos_plot";
    String LINX_FUSION_TSV = "linx_fusion_tsv";
    String LINX_BREAKEND_TSV = "linx_breakend_tsv";
    String LINX_SVS_TSV = "linx_svs_tsv";
    String LINX_DRIVER_CATALOG_TSV = "linx_driver_catalog_tsv";
    String CHORD_PREDICTION_TXT = "chord_prediction_txt";
    String MOLECULAR_TISSUE_ORIGIN_TXT = "molecular_tissue_origin_txt";
    String MOLECULAR_TISSUE_ORIGIN_PLOT = "molecular_tissue_origin_plot";
    String ANNOTATED_VIRUS_TSV = "annotated_virus_tsv";
    String PROTECT_EVIDENCE_TSV = "protect_evidence_tsv";

    // Resources used for generating an analysed patient report
    String GERMLINE_REPORTING_TSV = "germline_reporting_tsv";
    String SAMPLE_SUMMARY_TSV = "sample_summary_tsv";

    // Some additional optional params and flags
    String COMMENTS = "comments";
    String CORRECTED_REPORT = "corrected_report";
    String CORRECTED_REPORT_EXTERN = "corrected_report_extern";
    String LOG_DEBUG = "log_debug";
    String ONLY_CREATE_PDF = "only_create_pdf";

    // parameters for pipeline version
    String REQUIRE_PIPELINE_VERSION_FILE = "require_pipeline_version_file";
    String PIPELINE_VERSION_FILE = "pipeline_version_file";
    String EXPECTED_PIPELINE_VERSION = "expected_pipeline_version";
    String OVERRIDE_PIPELINE_VERSION = "override_pipeline_version";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which a patient report will be generated.");
        options.addOption(TUMOR_SAMPLE_BARCODE, true, "The sample barcode for which a patient report will be generated.");
        options.addOption(OUTPUT_DIRECTORY_REPORT, true, "Path to where the PDF report will be written to.");
        options.addOption(OUTPUT_DIRECTORY_DATA, true, "Path to where the data of the report will be written to.");

        options.addOption(PRIMARY_TUMOR_TSV, true, "Path towards the (curated) primary tumor TSV.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the directory holding the LIMS data");
        options.addOption(PEACH_GENOTYPE_TSV, true, "Path towards the peach genotype TSV.");

        options.addOption(RVA_LOGO, true, "Path towards an image file containing the RVA logo.");
        options.addOption(COMPANY_LOGO, true, "Path towards an image file containing the company logo.");
        options.addOption(SIGNATURE, true, "Path towards an image file containing the signature to be appended at the end of the report.");

        options.addOption(UDI_DI, true, "Code of the UDI DI code");

        options.addOption(REF_SAMPLE_ID, true, "The reference sample ID for the tumor sample for which a report is generated.");
        options.addOption(REF_SAMPLE_BARCODE, true, "The reference sample barcode for the tumor sample for which a report is generated.");

        options.addOption(QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(QC_FAIL_REASON, true, "One of: " + Strings.join(Lists.newArrayList(QCFailReason.validIdentifiers()), ','));

        options.addOption(PANEL, false, "If set, generates a panel report.");
        options.addOption(PANEL_QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(PANEL_QC_FAIL_REASON, true, "One of: " + Strings.join(Lists.newArrayList(QCFailReason.validIdentifiers()), ','));
        options.addOption(PANEL_VCF_NAME, true, "The name of the VCF file of the panel results.");

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");
        options.addOption(PURPLE_SOMATIC_DRIVER_CATALOG_TSV, true, "Path towards the purple somatic driver catalog TSV.");
        options.addOption(PURPLE_GERMLINE_DRIVER_CATALOG_TSV, true, "Path towards the purple germline driver catalog TSV.");
        options.addOption(PURPLE_SOMATIC_VARIANT_VCF, true, "Path towards the purple somatic variant VCF.");
        options.addOption(PURPLE_GERMLINE_VARIANT_VCF, true, "Path towards the purple germline variant VCF.");
        options.addOption(PURPLE_SOMATIC_COPYNUMBER_TSV, true, "Path towards the purple somatic copynumber TSV.");
        options.addOption(PURPLE_CIRCOS_PLOT, true, "Path towards the purple circos plot.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(LINX_BREAKEND_TSV, true, "Path towards the linx breakend TSV.");
        options.addOption(LINX_SVS_TSV, true, "Path towards the linx svs TSV.");
        options.addOption(LINX_DRIVER_CATALOG_TSV, true, "Path towards the LINX driver catalog TSV.");
        options.addOption(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT.");
        options.addOption(MOLECULAR_TISSUE_ORIGIN_TXT, true, "Path towards the molecular tissue origin TXT.");
        options.addOption(MOLECULAR_TISSUE_ORIGIN_PLOT, true, "Path towards the molecular tissue origin plot.");
        options.addOption(ANNOTATED_VIRUS_TSV, true, "Path towards the annotated virus TSV.");
        options.addOption(PROTECT_EVIDENCE_TSV, true, "Path towards the protect evidence TSV.");

        options.addOption(GERMLINE_REPORTING_TSV, true, "Path towards a TSV containing germline reporting config.");
        options.addOption(SAMPLE_SUMMARY_TSV, true, "Path towards a TSV containing the (clinical) summaries of the samples.");

        options.addOption(COMMENTS, true, "Additional comments to be added to the report (optional).");
        options.addOption(CORRECTED_REPORT, false, "If provided, generate a corrected report with corrected name");
        options.addOption(CORRECTED_REPORT_EXTERN, false, "If provided, generate a corrected report with intern/extern correction");

        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");
        options.addOption(ONLY_CREATE_PDF, false, "If provided, just the PDF will be generated and no additional data will be updated.");

        options.addOption(REQUIRE_PIPELINE_VERSION_FILE, false, "Boolean for determine pipeline version file is requierde");
        options.addOption(PIPELINE_VERSION_FILE, true, "Path towards the pipeline version (optional)");
        options.addOption(EXPECTED_PIPELINE_VERSION, true, "String of the expected pipeline version");
        options.addOption(OVERRIDE_PIPELINE_VERSION, false, "if set, the check for pipeline version is overridden");

        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version to use (either '37' or '38')");

        return options;
    }

    @Nullable
    String refSampleId();

    @Nullable
    String refSampleBarcode();

    @NotNull
    String tumorSampleId();

    @NotNull
    String tumorSampleBarcode();

    @NotNull
    String outputDirReport();

    @NotNull
    String outputDirData();

    @NotNull
    String primaryTumorTsv();

    @NotNull
    String limsDir();

    @NotNull
    String peachGenotypeTsv();

    @NotNull
    String rvaLogo();

    @NotNull
    String companyLogo();

    @NotNull
    String signature();

    @NotNull
    String udiDi();

    boolean qcFail();

    @Nullable
    QCFailReason qcFailReason();

    boolean panel();

    boolean panelQcFail();

    @NotNull
    String panelVCFname();

    @Nullable
    PanelFailReason panelQcFailReason();

    @NotNull
    String purplePurityTsv();

    @NotNull
    String purpleQcFile();

    @NotNull
    String purpleSomaticDriverCatalogTsv();

    @NotNull
    String purpleGermlineDriverCatalogTsv();

    @NotNull
    String purpleSomaticVariantVcf();

    @NotNull
    String purpleGermlineVariantVcf();

    @NotNull
    String purpleSomaticCopyNumberTsv();

    @NotNull
    String purpleCircosPlot();

    @NotNull
    String linxFusionTsv();

    @NotNull
    String linxBreakendTsv();

    @NotNull
    String linxSvsTsv();

    @NotNull
    String linxDriverCatalogTsv();

    @NotNull
    String chordPredictionTxt();

    @NotNull
    String molecularTissueOriginTxt();

    @NotNull
    String molecularTissueOriginPlot();

    @NotNull
    String annotatedVirusTsv();

    @NotNull
    String protectEvidenceTsv();

    @NotNull
    String germlineReportingTsv();

    @NotNull
    String sampleSummaryTsv();

    @Nullable
    String comments();

    boolean isCorrectedReport();

    boolean isCorrectedReportExtern();

    boolean onlyCreatePDF();

    boolean requirePipelineVersionFile();

    @Nullable
    String pipelineVersionFile();

    @NotNull
    String expectedPipelineVersion();

    boolean overridePipelineVersion();

    @NotNull
    RefGenomeVersion refGenomeVersion();

    @NotNull
    static PatientReporterConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
            LOGGER.debug("Switched root level logging to DEBUG");
        }

        boolean isQCFail = cmd.hasOption(QC_FAIL);
        boolean requirePipelineVersion = cmd.hasOption(REQUIRE_PIPELINE_VERSION_FILE);
        QCFailReason qcFailReason = null;
        if (isQCFail) {
            String qcFailReasonString = nonOptionalValue(cmd, QC_FAIL_REASON);
            qcFailReason = QCFailReason.fromIdentifier(qcFailReasonString);
            if (qcFailReason == null) {
                throw new ParseException("Did not recognize QC Fail reason: " + qcFailReasonString);
            }
        }

        String panelVCFFile = Strings.EMPTY;
        String pipelineVersion = null;

        boolean isPanel = cmd.hasOption(PANEL);
        boolean isPanelQCFail = cmd.hasOption(PANEL_QC_FAIL);
        PanelFailReason panelQcFailReason = null;
        if (isPanel) {
            if (isPanelQCFail) {
                String qcFailReasonString = nonOptionalValue(cmd, PANEL_QC_FAIL_REASON);
                panelQcFailReason = PanelFailReason.fromIdentifier(qcFailReasonString);
                if (panelQcFailReason == null) {
                    throw new ParseException("Did not recognize QC Fail reason: " + qcFailReasonString);
                }
            }
            if (requirePipelineVersion) {
                pipelineVersion = nonOptionalFile(cmd, PIPELINE_VERSION_FILE);
            }

            panelVCFFile= nonOptionalValue(cmd, PANEL_VCF_NAME);
        }


        String purplePurityTsv = Strings.EMPTY;
        String purpleQcFile = Strings.EMPTY;
        String purpleSomaticDriverCatalogTsv = Strings.EMPTY;
        String purpleGermlineDriverCatalogTsv = Strings.EMPTY;
        String purpleSomaticVariantVcf = Strings.EMPTY;
        String purpleGermlineVariantVcf = Strings.EMPTY;
        String purpleSomaticCopyNumberTsv = Strings.EMPTY;
        String purpleCircosPlot = Strings.EMPTY;
        String linxFusionTsv = Strings.EMPTY;
        String linxBreakendTsv = Strings.EMPTY;
        String linxSvsTsv = Strings.EMPTY;
        String linxDriverCatalogTsv = Strings.EMPTY;
        String chordPredictionTxt = Strings.EMPTY;
        String molecularTissueOriginTxt = Strings.EMPTY;
        String molecularTissueOriginPlot = Strings.EMPTY;
        String annotatedVirusTsv = Strings.EMPTY;
        String peachGenotypeTsv = Strings.EMPTY;
        String protectEvidenceTsv = Strings.EMPTY;

        String germlineReportingTsv = Strings.EMPTY;
        String sampleSummaryTsv = Strings.EMPTY;

        if (isQCFail && qcFailReason.isDeepWGSDataAvailable()) {
            if (requirePipelineVersion) {
                pipelineVersion = nonOptionalFile(cmd, PIPELINE_VERSION_FILE);
            }
            purplePurityTsv = nonOptionalFile(cmd, PURPLE_PURITY_TSV);
            purpleQcFile = nonOptionalFile(cmd, PURPLE_QC_FILE);
            peachGenotypeTsv = nonOptionalFile(cmd, PEACH_GENOTYPE_TSV);
        } else if (!isQCFail) {
            if (requirePipelineVersion) {
                pipelineVersion = nonOptionalFile(cmd, PIPELINE_VERSION_FILE);
            }

            purplePurityTsv = nonOptionalFile(cmd, PURPLE_PURITY_TSV);
            purpleQcFile = nonOptionalFile(cmd, PURPLE_QC_FILE);
            purpleSomaticDriverCatalogTsv = nonOptionalFile(cmd, PURPLE_SOMATIC_DRIVER_CATALOG_TSV);
            purpleGermlineDriverCatalogTsv = nonOptionalFile(cmd, PURPLE_GERMLINE_DRIVER_CATALOG_TSV);
            purpleSomaticVariantVcf = nonOptionalFile(cmd, PURPLE_SOMATIC_VARIANT_VCF);
            purpleGermlineVariantVcf = nonOptionalFile(cmd, PURPLE_GERMLINE_VARIANT_VCF);
            purpleSomaticCopyNumberTsv = nonOptionalFile(cmd, PURPLE_SOMATIC_COPYNUMBER_TSV);
            purpleCircosPlot = nonOptionalFile(cmd, PURPLE_CIRCOS_PLOT);
            linxFusionTsv = nonOptionalFile(cmd, LINX_FUSION_TSV);
            linxBreakendTsv = nonOptionalFile(cmd, LINX_BREAKEND_TSV);
            linxSvsTsv = nonOptionalFile(cmd, LINX_SVS_TSV);
            linxDriverCatalogTsv = nonOptionalFile(cmd, LINX_DRIVER_CATALOG_TSV);
            chordPredictionTxt = nonOptionalFile(cmd, CHORD_PREDICTION_TXT);
            molecularTissueOriginTxt = nonOptionalFile(cmd, MOLECULAR_TISSUE_ORIGIN_TXT);
            molecularTissueOriginPlot = nonOptionalFile(cmd, MOLECULAR_TISSUE_ORIGIN_PLOT);
            annotatedVirusTsv = nonOptionalFile(cmd, ANNOTATED_VIRUS_TSV);
            peachGenotypeTsv = nonOptionalFile(cmd, PEACH_GENOTYPE_TSV);
            protectEvidenceTsv = nonOptionalFile(cmd, PROTECT_EVIDENCE_TSV);

            germlineReportingTsv = nonOptionalFile(cmd, GERMLINE_REPORTING_TSV);
            sampleSummaryTsv = nonOptionalFile(cmd, SAMPLE_SUMMARY_TSV);
        }

        return ImmutablePatientReporterConfig.builder()
                .refSampleId(cmd.hasOption(REF_SAMPLE_ID) ? nonOptionalValue(cmd, REF_SAMPLE_ID) : null)
                .refSampleBarcode(cmd.hasOption(REF_SAMPLE_BARCODE) ? nonOptionalValue(cmd, REF_SAMPLE_BARCODE) : null)
                .tumorSampleId(nonOptionalValue(cmd, TUMOR_SAMPLE_ID))
                .tumorSampleBarcode(nonOptionalValue(cmd, TUMOR_SAMPLE_BARCODE))
                .outputDirReport(nonOptionalDir(cmd, OUTPUT_DIRECTORY_REPORT))
                .outputDirData(nonOptionalDir(cmd, OUTPUT_DIRECTORY_DATA))
                .primaryTumorTsv(nonOptionalFile(cmd, PRIMARY_TUMOR_TSV))
                .limsDir(nonOptionalDir(cmd, LIMS_DIRECTORY))
                .rvaLogo(nonOptionalFile(cmd, RVA_LOGO))
                .companyLogo(nonOptionalFile(cmd, COMPANY_LOGO))
                .signature(nonOptionalFile(cmd, SIGNATURE))
                .udiDi(nonOptionalValue(cmd, UDI_DI))
                .qcFail(isQCFail)
                .qcFailReason(qcFailReason)
                .panel(isPanel)
                .panelQcFail(isPanelQCFail)
                .panelQcFailReason(panelQcFailReason)
                .panelVCFname(panelVCFFile)
                .purplePurityTsv(purplePurityTsv)
                .purpleQcFile(purpleQcFile)
                .purpleSomaticDriverCatalogTsv(purpleSomaticDriverCatalogTsv)
                .purpleGermlineDriverCatalogTsv(purpleGermlineDriverCatalogTsv)
                .purpleSomaticVariantVcf(purpleSomaticVariantVcf)
                .purpleGermlineVariantVcf(purpleGermlineVariantVcf)
                .purpleSomaticCopyNumberTsv(purpleSomaticCopyNumberTsv)
                .purpleCircosPlot(purpleCircosPlot)
                .linxFusionTsv(linxFusionTsv)
                .linxBreakendTsv(linxBreakendTsv)
                .linxSvsTsv(linxSvsTsv)
                .linxDriverCatalogTsv(linxDriverCatalogTsv)
                .chordPredictionTxt(chordPredictionTxt)
                .molecularTissueOriginTxt(molecularTissueOriginTxt)
                .molecularTissueOriginPlot(molecularTissueOriginPlot)
                .annotatedVirusTsv(annotatedVirusTsv)
                .peachGenotypeTsv(peachGenotypeTsv)
                .protectEvidenceTsv(protectEvidenceTsv)
                .germlineReportingTsv(germlineReportingTsv)
                .sampleSummaryTsv(sampleSummaryTsv)
                .comments(cmd.getOptionValue(COMMENTS))
                .isCorrectedReport(cmd.hasOption(CORRECTED_REPORT))
                .isCorrectedReportExtern(cmd.hasOption(CORRECTED_REPORT_EXTERN))
                .onlyCreatePDF(cmd.hasOption(ONLY_CREATE_PDF))
                .requirePipelineVersionFile(requirePipelineVersion)
                .pipelineVersionFile(pipelineVersion)
                .expectedPipelineVersion(cmd.getOptionValue(EXPECTED_PIPELINE_VERSION))
                .overridePipelineVersion(cmd.hasOption(OVERRIDE_PIPELINE_VERSION))
                .refGenomeVersion(RefGenomeVersion.from(nonOptionalValue(cmd, RefGenomeVersion.REF_GENOME_VERSION)))
                .build();
    }

    @NotNull
    static String nonOptionalValue(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            throw new ParseException("Parameter must be provided: " + param);
        }

        return value;
    }

    @NotNull
    static String nonOptionalDir(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = nonOptionalValue(cmd, param);

        if (!pathExists(value) || !pathIsDirectory(value)) {
            throw new ParseException("Parameter '" + param + "' must be an existing directory: " + value);
        }

        return value;
    }

    @NotNull
    static String nonOptionalFile(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = nonOptionalValue(cmd, param);

        if (!pathExists(value)) {
            throw new ParseException("Parameter '" + param + "' must be an existing file: " + value);
        }

        return value;
    }

    static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }
}