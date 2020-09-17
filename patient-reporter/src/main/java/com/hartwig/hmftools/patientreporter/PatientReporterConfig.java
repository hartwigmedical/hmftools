package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.common.cli.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.cli.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;

import java.io.File;
import java.nio.file.Files;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PatientReporterConfig {

    // General params needed for every report
    String TUMOR_SAMPLE_ID = "tumor_sample_id";
    String TUMOR_SAMPLE_BARCODE = "tumor_sample_barcode";
    String OUTPUT_DIRECTORY_REPORT = "output_dir_report";
    String OUTPUT_DIRECTORY_DATA_REPORT = "output_dir_data_report";

    String REPORTING_DB_TSV = "reporting_db_tsv";
    String TUMOR_LOCATION_CSV = "tumor_location_csv";
    String LIMS_DIRECTORY = "lims_dir";

    String RVA_LOGO = "rva_logo";
    String COMPANY_LOGO = "company_logo";
    String SIGNATURE = "signature";

    // General params needed for every report but for QC fail it could be null
    String REF_SAMPLE_ID = "ref_sample_id";
    String REF_SAMPLE_BARCODE = "ref_sample_barcode";

    // Params specific for QC Fail reports
    String QC_FAIL = "qc_fail";
    String QC_FAIL_REASON = "qc_fail_reason";

    // Params specific for actual patient reports
    String PURPLE_PURITY_TSV = "purple_purity_tsv";
    String PURPLE_QC_FILE = "purple_qc_file";
    String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    String BACHELOR_TSV = "bachelor_tsv";
    String LINX_FUSION_TSV = "linx_fusion_tsv";
    String LINX_DISRUPTION_TSV = "linx_disruption_tsv";
    String LINX_VIRAL_INSERTION_TSV = "linx_viral_insertion_tsv";
    String LINX_DRIVERS_TSV = "linx_drivers_tsv";
    String CHORD_PREDICTION_TXT = "chord_prediction_txt";
    String CIRCOS_FILE = "circos_file";

    String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    String GERMLINE_GENES_CSV = "germline_genes_csv";
    String SAMPLE_SUMMARY_TSV = "sample_summary_tsv";

    // Some additional optional params and flags
    String COMMENTS = "comments";
    String CORRECTED_REPORT = "corrected_report";
    String UNOFFICIAL_REPORT = "unofficial_report";
    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);

        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which a patient report will be generated.");
        options.addOption(TUMOR_SAMPLE_BARCODE, true, "The sample barcode for which a patient report will be generated.");
        options.addOption(OUTPUT_DIRECTORY_REPORT, true, "Path to where the PDF report will be written to.");
        options.addOption(OUTPUT_DIRECTORY_DATA_REPORT, true, "Path to where the data of the report will be written to.");

        options.addOption(REPORTING_DB_TSV, true, "Path towards output file for the reporting db TSV.");
        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the directory holding the LIMS data");

        options.addOption(RVA_LOGO, true, "Path towards an image file containing the RVA logo.");
        options.addOption(COMPANY_LOGO, true, "Path towards an image file containing the company logo.");
        options.addOption(SIGNATURE, true, "Path towards an image file containing the signature to be appended at the end of the report.");

        options.addOption(REF_SAMPLE_ID, false, "The reference sample ID for the tumor sample for which we are generating a report.");
        options.addOption(REF_SAMPLE_BARCODE,
                false,
                "The reference sample barcode for the tumor sample for which we are generating a report.");

        options.addOption(QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(QC_FAIL_REASON, true, "One of: " + Strings.join(Lists.newArrayList(QCFailReason.validIdentifiers()), ','));

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(BACHELOR_TSV, true, "Path towards the germline TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(LINX_DISRUPTION_TSV, true, "Path towards the linx disruption TSV.");
        options.addOption(LINX_VIRAL_INSERTION_TSV, true, "Path towards the LINX viral insertion TSV.");
        options.addOption(LINX_DRIVERS_TSV, true, "Path towards the LINX driver catalog TSV.");
        options.addOption(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT.");
        options.addOption(CIRCOS_FILE, true, "Path towards the circos file.");

        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the directory holding knowledgebase output files.");
        options.addOption(GERMLINE_GENES_CSV, true, "Path towards a CSV containing germline genes which we want to report.");
        options.addOption(SAMPLE_SUMMARY_TSV, true, "Path towards a TSV containing the (clinical) summaries of the samples.");

        options.addOption(COMMENTS, true, "Additional comments to be added to the report (optional).");
        options.addOption(CORRECTED_REPORT, false, "If provided, generate a corrected report with corrected name");
        options.addOption(UNOFFICIAL_REPORT, false, "If provided, generates a report with potentially some sections altered.");
        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");
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
    String reportingDbTsv();

    @NotNull
    String tumorLocationCsv();

    @NotNull
    String limsDir();

    @NotNull
    String rvaLogo();

    @NotNull
    String companyLogo();

    @NotNull
    String signature();

    boolean qcFail();

    @NotNull
    QCFailReason qcFailReason();

    @NotNull
    String purplePurityTsv();

    @NotNull
    String purpleQcFile();

    @NotNull
    String purpleGeneCnvTsv();

    @NotNull
    String somaticVariantVcf();

    @NotNull
    String bachelorTsv();

    @NotNull
    String linxFusionTsv();

    @NotNull
    String linxDisruptionTsv();

    @NotNull
    String linxViralInsertionTsv();

    @NotNull
    String linxDriversTsv();

    @NotNull
    String chordPredictionTxt();

    @NotNull
    String circosFile();

    @NotNull
    String knowledgebaseDir();

    @NotNull
    String germlineGenesCsv();

    @NotNull
    String sampleSummaryTsv();

    @NotNull
    String driverGenePanelTsv();

    @Nullable
    String comments();

    boolean correctedReport();

    boolean unofficialReport();

    @NotNull
    static PatientReporterConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        boolean isQCFail = cmd.hasOption(QC_FAIL);
        QCFailReason qcFailReason = QCFailReason.UNDEFINED;
        if (isQCFail) {
            String qcFailReasonString = nonOptionalValue(cmd, QC_FAIL_REASON);
            qcFailReason = QCFailReason.fromIdentifier(qcFailReasonString);
            if (qcFailReason == QCFailReason.UNDEFINED) {
                throw new ParseException("Did not recognize QC Fail reason: " + qcFailReasonString);
            }
        }

        String purplePurityTsv = Strings.EMPTY;
        String purpleQCFile = Strings.EMPTY;
        String purpleGeneCnvTsv = Strings.EMPTY;
        String somaticVariantVcf = Strings.EMPTY;
        String bachelorTsv = Strings.EMPTY;
        String linxFusionTsv = Strings.EMPTY;
        String linxDisruptionTsv = Strings.EMPTY;
        String linxViralInsertionTsv = Strings.EMPTY;
        String linxDriversTsv = Strings.EMPTY;
        String chordPredictionTxt = Strings.EMPTY;
        String circosFile = Strings.EMPTY;
        String knowledgebaseDirectory = Strings.EMPTY;
        String germlineGenesCsv = Strings.EMPTY;
        String sampleSummaryTsv = Strings.EMPTY;

        if (qcFailReason.isDeepWGSDataAvailable()) {
            purplePurityTsv = nonOptionalFile(cmd, PURPLE_PURITY_TSV);
        } else if (!isQCFail) {
            purpleQCFile = nonOptionalFile(cmd, PURPLE_QC_FILE);
            purplePurityTsv = nonOptionalFile(cmd, PURPLE_PURITY_TSV);
            purpleGeneCnvTsv = nonOptionalFile(cmd, PURPLE_GENE_CNV_TSV);
            somaticVariantVcf = nonOptionalFile(cmd, SOMATIC_VARIANT_VCF);
            bachelorTsv = nonOptionalFile(cmd, BACHELOR_TSV);
            linxFusionTsv = nonOptionalFile(cmd, LINX_FUSION_TSV);
            linxDisruptionTsv = nonOptionalFile(cmd, LINX_DISRUPTION_TSV);
            linxViralInsertionTsv = nonOptionalFile(cmd, LINX_VIRAL_INSERTION_TSV);
            linxDriversTsv = nonOptionalFile(cmd, LINX_DRIVERS_TSV);
            chordPredictionTxt = nonOptionalFile(cmd, CHORD_PREDICTION_TXT);
            circosFile = nonOptionalFile(cmd, CIRCOS_FILE);
            knowledgebaseDirectory = nonOptionalDir(cmd, KNOWLEDGEBASE_DIRECTORY);
            germlineGenesCsv = nonOptionalFile(cmd, GERMLINE_GENES_CSV);
            sampleSummaryTsv = nonOptionalFile(cmd, SAMPLE_SUMMARY_TSV);
        }

        return ImmutablePatientReporterConfig.builder()
                .refSampleId(cmd.hasOption(REF_SAMPLE_ID) ? optionalValue(cmd, REF_SAMPLE_ID) : null )
                .refSampleBarcode(cmd.hasOption(REF_SAMPLE_BARCODE) ? optionalValue(cmd, REF_SAMPLE_BARCODE): null)
                .tumorSampleId(nonOptionalValue(cmd, TUMOR_SAMPLE_ID))
                .tumorSampleBarcode(nonOptionalValue(cmd, TUMOR_SAMPLE_BARCODE))
                .outputDirReport(nonOptionalDir(cmd, OUTPUT_DIRECTORY_REPORT))
                .outputDirData(nonOptionalDir(cmd, OUTPUT_DIRECTORY_DATA_REPORT))
                .reportingDbTsv(nonOptionalFile(cmd, REPORTING_DB_TSV))
                .tumorLocationCsv(nonOptionalFile(cmd, TUMOR_LOCATION_CSV))
                .limsDir(nonOptionalDir(cmd, LIMS_DIRECTORY))
                .rvaLogo(nonOptionalFile(cmd, RVA_LOGO))
                .companyLogo(nonOptionalFile(cmd, COMPANY_LOGO))
                .signature(nonOptionalFile(cmd, SIGNATURE))
                .driverGenePanelTsv(nonOptionalFile(cmd, DRIVER_GENE_PANEL_OPTION))
                .qcFail(isQCFail)
                .qcFailReason(qcFailReason)
                .purplePurityTsv(purplePurityTsv)
                .purpleQcFile(purpleQCFile)
                .purpleGeneCnvTsv(purpleGeneCnvTsv)
                .somaticVariantVcf(somaticVariantVcf)
                .bachelorTsv(bachelorTsv)
                .linxFusionTsv(linxFusionTsv)
                .linxDisruptionTsv(linxDisruptionTsv)
                .linxViralInsertionTsv(linxViralInsertionTsv)
                .linxDriversTsv(linxDriversTsv)
                .chordPredictionTxt(chordPredictionTxt)
                .circosFile(circosFile)
                .knowledgebaseDir(knowledgebaseDirectory)
                .germlineGenesCsv(germlineGenesCsv)
                .sampleSummaryTsv(sampleSummaryTsv)
                .comments(cmd.getOptionValue(COMMENTS))
                .correctedReport(cmd.hasOption(CORRECTED_REPORT))
                .unofficialReport(cmd.hasOption(UNOFFICIAL_REPORT))
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

    @Nullable
    static String optionalValue(@NotNull CommandLine cmd, @NotNull String param) {
        return cmd.getOptionValue(param);
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
