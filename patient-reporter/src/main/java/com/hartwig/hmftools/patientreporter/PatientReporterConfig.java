package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.nio.file.Files;

import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface PatientReporterConfig {
    Logger LOGGER = LogManager.getLogger(PatientReporterConfig.class);

    // General params needed for every report
    String REF_SAMPLE_ID = "ref_sample_id";
    String REF_SAMPLE_BARCODE = "ref_sample_barcode";
    String TUMOR_SAMPLE_ID = "tumor_sample_id";
    String TUMOR_SAMPLE_BARCODE = "tumor_sample_barcode";
    String OUTPUT_DIRECTORY = "output_dir";

    String REPORTING_DB_TSV = "reporting_db_tsv";
    String TUMOR_LOCATION_CSV = "tumor_location_csv";
    String LIMS_DIRECTORY = "lims_dir";
    String HOSPITAL_DIRECTORY = "hospital_dir";
    String CONTACT_WIDE_TSV = "contact_wide_tsv";

    String RVA_LOGO = "rva_logo";
    String COMPANY_LOGO = "company_logo";
    String SIGNATURE = "signature";

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
        options.addOption(REF_SAMPLE_ID, true, "The reference sample ID for the sample for which we are generating a report.");
        options.addOption(REF_SAMPLE_BARCODE, true, "The reference sample barcode for the sample for which we are generating a report.");
        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which a patient report will be generated.");
        options.addOption(TUMOR_SAMPLE_BARCODE, true, "The sample barcode for which a patient report will be generated.");
        options.addOption(OUTPUT_DIRECTORY, true, "Path to where the PDF report will be written to.");

        options.addOption(REPORTING_DB_TSV, true, "Path towards output file for the reporting db TSV.");
        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the directory holding the LIMS data");
        options.addOption(HOSPITAL_DIRECTORY, true, "Path towards the directory containing hospital data.");
        options.addOption(CONTACT_WIDE_TSV, true, "Path towards the file of contact for WIDE TSV.");

        options.addOption(RVA_LOGO, true, "Path towards a image file containing the RVA logo.");
        options.addOption(COMPANY_LOGO, true, "Path towards a image file containing the company logo.");
        options.addOption(SIGNATURE, true, "Path towards a image file containing the signature to be appended at the end of the report.");

        options.addOption(QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(QC_FAIL_REASON,
                true,
                "Either 'low_dna_yield', 'post_analysis_fail', 'shallow_seq' or 'insufficient_tissue_delivered'");

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(BACHELOR_TSV, true, "Path towards the germline TSV (optional).");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(LINX_DISRUPTION_TSV, true, "Path towards the linx disruption TSV.");
        options.addOption(LINX_VIRAL_INSERTION_TSV, true, "Path towards the LINX viral insertion TSV.");
        options.addOption(LINX_DRIVERS_TSV, true, "Path towards the LINX driver catalog TSV.");
        options.addOption(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT .");
        options.addOption(CIRCOS_FILE, true, "Path towards the circos file.");

        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the directory holding knowledgebase output files.");
        options.addOption(GERMLINE_GENES_CSV, true, "Path towards a CSV containing germline genes which we want to report.");
        options.addOption(SAMPLE_SUMMARY_TSV, true, "Path towards a TSV containing the (clinical) summaries of the samples.");

        options.addOption(COMMENTS, true, "Additional comments to be added to the report (optional).");
        options.addOption(CORRECTED_REPORT, false, "If provided, generate a corrected report with corrected name");
        options.addOption(UNOFFICIAL_REPORT, false, "If provided, generates a report with potentially some sections removed.");
        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");
        return options;
    }

    @NotNull
    String refSampleID();

    @NotNull
    String refSampleBarcode();

    @NotNull
    String tumorSampleId();

    @NotNull
    String tumorSampleBarcode();

    @NotNull
    String outputDir();

    @NotNull
    String reportingDbTsv();

    @NotNull
    String tumorLocationCsv();

    @NotNull
    String limsDir();

    @NotNull
    String hospitalDir();

    @NotNull
    String contactWideTsv();

    @NotNull
    String rvaLogo();

    @NotNull
    String companyLogo();

    @NotNull
    String signature();

    @NotNull
    String qcFail();

    @NotNull
    String qcFailReason();

    @NotNull
    String purplePurityTsv();

    @NotNull
    String purpleQcFile();

    @NotNull
    String PurpleGeneCnvTsv();

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
    String comments();

    @NotNull
    String correctedReport();

    @NotNull
    String unofficialReport();

    @NotNull
    String logDebug();

    @NotNull
    static PatientReporterConfig createConfig(@NotNull final String version, @NotNull final CommandLine cmd) throws ParseException {
        return ImmutablePatientReporterConfig.builder()
                .refSampleID(cmd.getOptionValue(REF_SAMPLE_ID))
                .refSampleBarcode(cmd.getOptionValue(REF_SAMPLE_BARCODE))
                .tumorSampleId(cmd.getOptionValue(TUMOR_SAMPLE_ID))
                .tumorSampleBarcode(cmd.getOptionValue(TUMOR_SAMPLE_BARCODE))
                .outputDir(cmd.getOptionValue(OUTPUT_DIRECTORY))
                .reportingDbTsv(cmd.getOptionValue(REPORTING_DB_TSV))
                .tumorLocationCsv(cmd.getOptionValue(TUMOR_LOCATION_CSV))
                .limsDir(cmd.getOptionValue(LIMS_DIRECTORY))
                .hospitalDir(cmd.getOptionValue(HOSPITAL_DIRECTORY))
                .contactWideTsv(cmd.getOptionValue(CONTACT_WIDE_TSV))
                .rvaLogo(cmd.getOptionValue(RVA_LOGO))
                .companyLogo(cmd.getOptionValue(COMPANY_LOGO))
                .signature(cmd.getOptionValue(SIGNATURE))
                .qcFail(cmd.getOptionValue(QC_FAIL))
                .qcFailReason(cmd.getOptionValue(QC_FAIL_REASON))
                .purplePurityTsv(cmd.getOptionValue(PURPLE_PURITY_TSV))
                .purpleQcFile(cmd.getOptionValue(PURPLE_QC_FILE))
                .PurpleGeneCnvTsv(cmd.getOptionValue(PURPLE_GENE_CNV_TSV))
                .somaticVariantVcf(cmd.getOptionValue(SOMATIC_VARIANT_VCF))
                .bachelorTsv(cmd.getOptionValue(BACHELOR_TSV))
                .linxFusionTsv(cmd.getOptionValue(LINX_FUSION_TSV))
                .linxDisruptionTsv(cmd.getOptionValue(LINX_DISRUPTION_TSV))
                .linxViralInsertionTsv(cmd.getOptionValue(LINX_VIRAL_INSERTION_TSV))
                .linxDriversTsv(cmd.getOptionValue(LINX_DRIVERS_TSV))
                .chordPredictionTxt(cmd.getOptionValue(CHORD_PREDICTION_TXT))
                .circosFile(cmd.getOptionValue(CIRCOS_FILE))
                .knowledgebaseDir(cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY))
                .germlineGenesCsv(cmd.getOptionValue(GERMLINE_GENES_CSV))
                .sampleSummaryTsv(cmd.getOptionValue(SAMPLE_SUMMARY_TSV))
                .comments(cmd.getOptionValue(COMMENTS))
                .correctedReport(cmd.getOptionValue(CORRECTED_REPORT))
                .unofficialReport(cmd.getOptionValue(UNOFFICIAL_REPORT))
                .logDebug(cmd.getOptionValue(LOG_DEBUG))
                .build();

    }

    static boolean validInputForBaseReport(@NotNull CommandLine cmd) {
        return valueExists(cmd, REF_SAMPLE_ID) && valueExists(cmd, REF_SAMPLE_BARCODE) && valueExists(cmd, TUMOR_SAMPLE_ID) && valueExists(
                cmd,
                TUMOR_SAMPLE_BARCODE) && dirExists(cmd, OUTPUT_DIRECTORY) && fileExists(cmd, REPORTING_DB_TSV) && fileExists(cmd,
                TUMOR_LOCATION_CSV) && dirExists(cmd, LIMS_DIRECTORY) && dirExists(cmd, HOSPITAL_DIRECTORY) && fileExists(cmd, SIGNATURE)
                && fileExists(cmd, RVA_LOGO) && fileExists(cmd, COMPANY_LOGO) && fileExists(cmd, CONTACT_WIDE_TSV);
    }

    static boolean validInputForAnalysedReport(@NotNull CommandLine cmd) {
        return fileExists(cmd, PURPLE_PURITY_TSV) && fileExists(cmd, PURPLE_QC_FILE) && fileExists(cmd, PURPLE_GENE_CNV_TSV) && fileExists(
                cmd,
                SOMATIC_VARIANT_VCF) && fileExists(cmd, BACHELOR_TSV) && fileExists(cmd, LINX_FUSION_TSV) && fileExists(cmd,
                LINX_DISRUPTION_TSV) && fileExists(cmd, LINX_VIRAL_INSERTION_TSV) && fileExists(cmd, LINX_DRIVERS_TSV) && fileExists(cmd,
                CHORD_PREDICTION_TXT) && fileExists(cmd, CIRCOS_FILE) && dirExists(cmd, KNOWLEDGEBASE_DIRECTORY) && fileExists(cmd,
                GERMLINE_GENES_CSV) && fileExists(cmd, SAMPLE_SUMMARY_TSV);
    }

    static boolean validInputForQCFailReport(@NotNull CommandLine cmd) {
        QCFailReason qcFailReason = QCFailReason.fromIdentifier(cmd.getOptionValue(QC_FAIL_REASON));
        if (qcFailReason == QCFailReason.UNDEFINED) {
            LOGGER.warn("{} has to be 'low_dna_yield', 'post_analysis_fail', 'shallow_seq_low_purity' or 'insufficient_tissue_delivered'",
                    QC_FAIL_REASON);
        } else {
            return true;
        }
        return false;
    }

    static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            LOGGER.warn("'{}' has to be provided", param);
            return false;
        }
        return true;
    }

    static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value)) {
            LOGGER.warn("'{}' has to be an existing file: '{}'", param, value);
            return false;
        }

        return true;
    }

    static boolean dirExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value) || !pathIsDirectory(value)) {
            LOGGER.warn("'{}' has to be an existing directory: '{}'", param, value);
            return false;
        }

        return true;
    }

    static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }

}
