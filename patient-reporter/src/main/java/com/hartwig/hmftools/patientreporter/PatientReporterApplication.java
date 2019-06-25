package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.hospital.HospitalModel;
import com.hartwig.hmftools.common.hospital.HospitalModelFactory;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientreporter.cfreport.CFReportWriter;
import com.hartwig.hmftools.patientreporter.qcfail.ImmutableQCFailReporter;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReporter;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);

    public static final String VERSION = PatientReporterApplication.class.getPackage().getImplementationVersion();

    // Uncomment this line when generating an example report using PDFWriterTest
    //        public static final String VERSION = "7.0";

    // General params needed for every report
    private static final String TUMOR_SAMPLE = "tumor_sample";
    private static final String OUTPUT_DIRECTORY = "output_dir";

    private static final String LIMS_DIRECTORY = "lims_dir";
    private static final String HOSPITAL_DIRECTORY = "hospital_dir";
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";

    private static final String RVA_LOGO = "rva_logo";
    private static final String COMPANY_LOGO = "company_logo";
    private static final String SIGNATURE = "signature";

    // Params specific for QC Fail reports
    private static final String QC_FAIL = "qc_fail";
    private static final String QC_FAIL_REASON = "qc_fail_reason";

    // Params specific for actual patient reports
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";
    private static final String LINX_DISRUPTION_TSV = "linx_disruption_tsv";
    private static final String BACHELOR_CSV = "bachelor_csv";
    private static final String CHORD_PREDICTION_FILE = "chord_prediction_file";
    private static final String CIRCOS_FILE = "circos_file";

    private static final String REF_SAMPLE = "ref_sample";
    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String DRUP_GENES_CSV = "drup_genes_csv";
    private static final String HOTSPOT_TSV = "hotspot_tsv";
    private static final String GERMLINE_GENES_CSV = "germline_genes_csv";
    private static final String SAMPLE_SUMMARY_CSV = "sample_summary_csv";
    private static final String FASTA_FILE_LOCATION = "fasta_file_location";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    // Some additional optional params
    private static final String COMMENTS = "comments";
    private static final String LOG_DEBUG = "log_debug";

    public static void main(final String... args) throws ParseException, IOException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        if (!validInputForReportWriter(cmd) || !validInputForBaseReportData(cmd)) {
            printUsageAndExit(options);
        }

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        LOGGER.info("Running patient reporter v" + VERSION);
        final ReportWriter reportWriter = CFReportWriter.createProductionReportWriter();

        if (cmd.hasOption(QC_FAIL) && validInputForQCFailReport(cmd)) {
            final String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
            LOGGER.info("Generating qc-fail report for {}", tumorSample);
            final QCFailReason reason = QCFailReason.fromIdentifier(cmd.getOptionValue(QC_FAIL_REASON));
            final QCFailReporter reporter = ImmutableQCFailReporter.of(buildBaseReportData(cmd));

            final QCFailReport report = reporter.run(tumorSample, reason, cmd.getOptionValue(COMMENTS));
            final String outputFilePath = generateOutputFilePathForPatientReport(cmd.getOptionValue(OUTPUT_DIRECTORY), report);
            reportWriter.writeQCFailReport(report, outputFilePath);
        } else if (validInputForAnalysedSample(cmd)) {
            final String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
            LOGGER.info("Generating patient report for {}", tumorSample);
            final SequencedReportData reporterData = buildReporterData(cmd);
            final PatientReporter reporter = buildReporter(cmd, reporterData);

            final AnalysedPatientReport report = reporter.run(tumorSample,
                    cmd.getOptionValue(REF_SAMPLE),
                    cmd.getOptionValue(PURPLE_PURITY_TSV),
                    cmd.getOptionValue(PURPLE_GENE_CNV_TSV),
                    cmd.getOptionValue(SOMATIC_VARIANT_VCF),
                    cmd.getOptionValue(BACHELOR_CSV),
                    cmd.getOptionValue(CHORD_PREDICTION_FILE),
                    cmd.getOptionValue(CIRCOS_FILE),
                    cmd.getOptionValue(COMMENTS));
            final String outputFilePath = generateOutputFilePathForPatientReport(cmd.getOptionValue(OUTPUT_DIRECTORY), report);
            reportWriter.writeAnalysedPatientReport(report, outputFilePath);
        } else {
            printUsageAndExit(options);
        }
    }

    @NotNull
    private static String generateOutputFilePathForPatientReport(@NotNull String reportDirectory, @NotNull PatientReport patientReport) {
        SampleReport sampleReport = patientReport.sampleReport();
        LimsSampleType type = LimsSampleType.fromSampleId(sampleReport.sampleId());

        String filePrefix =
                type == LimsSampleType.CORE ? sampleReport.sampleId() + "_" + sampleReport.hospitalPatientId() : sampleReport.sampleId();
        return reportDirectory + File.separator + filePrefix + "_hmf_report.pdf";
    }

    @NotNull
    private static BaseReportData buildBaseReportData(@NotNull final CommandLine cmd) throws IOException {
        String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        LOGGER.info("Loading tumor location CSV from {}.", tumorLocationCsv);
        final List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info(" Loaded tumor locations for {} patients.", patientTumorLocations.size());

        String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        LOGGER.info("Loading LIMS database from {}.", limsDirectory);
        final Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
        LOGGER.info(" Loaded data for {} samples.", lims.sampleCount());

        String hospitalsDirectory = cmd.getOptionValue(HOSPITAL_DIRECTORY);
        LOGGER.info("Loading hospitals from {}.", hospitalsDirectory);
        final HospitalModel hospitalModel = HospitalModelFactory.fromHospitalDirectory(hospitalsDirectory);
        LOGGER.info(" Loaded data for {} hospitals.", hospitalModel.hospitalCount());

        return ImmutableBaseReportData.of(patientTumorLocations,
                lims,
                hospitalModel,
                cmd.getOptionValue(SIGNATURE),
                cmd.getOptionValue(RVA_LOGO),
                cmd.getOptionValue(COMPANY_LOGO));
    }

    @NotNull
    private static SequencedReportData buildReporterData(@NotNull final CommandLine cmd) throws IOException {
        return SequencedReportDataLoader.buildFromFiles(cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY),
                cmd.getOptionValue(DRUP_GENES_CSV),
                cmd.getOptionValue(HOTSPOT_TSV),
                cmd.getOptionValue(FASTA_FILE_LOCATION),
                cmd.getOptionValue(HIGH_CONFIDENCE_BED),
                cmd.getOptionValue(GERMLINE_GENES_CSV),
                cmd.getOptionValue(SAMPLE_SUMMARY_CSV));
    }

    @NotNull
    private static PatientReporter buildReporter(@NotNull final CommandLine cmd, @NotNull final SequencedReportData sequencedReportData)
            throws IOException {
        final SvAnalyzer svAnalyzer = SvAnalyzer.fromFiles(cmd.getOptionValue(LINX_FUSION_TSV), cmd.getOptionValue(LINX_DISRUPTION_TSV));

        return new PatientReporter(buildBaseReportData(cmd), sequencedReportData, svAnalyzer);
    }

    private static boolean validInputForAnalysedSample(@NotNull final CommandLine cmd) {
        final String purplePurityTsv = cmd.getOptionValue(PURPLE_PURITY_TSV);
        final String purpleGeneCnvTsv = cmd.getOptionValue(PURPLE_GENE_CNV_TSV);
        final String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        final String linxFusionsTsv = cmd.getOptionValue(LINX_FUSION_TSV);
        final String linxDisruptionsTsv = cmd.getOptionValue(LINX_DISRUPTION_TSV);
        final String bachelorCsv = cmd.getOptionValue(BACHELOR_CSV);
        final String chordPredictionFile = cmd.getOptionValue(CHORD_PREDICTION_FILE);
        final String circosFile = cmd.getOptionValue(CIRCOS_FILE);

        final String refSample = cmd.getOptionValue(REF_SAMPLE);
        final String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        final String drupGenesCsv = cmd.getOptionValue(DRUP_GENES_CSV);
        final String hotspotTsv = cmd.getOptionValue(HOTSPOT_TSV);
        final String germlineGenesCsv = cmd.getOptionValue(GERMLINE_GENES_CSV);
        final String sampleSummaryCsv = cmd.getOptionValue(SAMPLE_SUMMARY_CSV);
        final String fastaFileLocation = cmd.getOptionValue(FASTA_FILE_LOCATION);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);

        if (purplePurityTsv == null || !exists(purplePurityTsv)) {
            LOGGER.warn(PURPLE_PURITY_TSV + " has to be an existing file: " + purplePurityTsv);
        } else if (purpleGeneCnvTsv == null || !exists(purpleGeneCnvTsv)) {
            LOGGER.warn(PURPLE_GENE_CNV_TSV + " has to be an existing file: " + purpleGeneCnvTsv);
        } else if (somaticVariantVcf == null || !exists(somaticVariantVcf)) {
            LOGGER.warn(SOMATIC_VARIANT_VCF + " has to be an existing file: " + somaticVariantVcf);
        } else if (linxFusionsTsv == null || !exists(linxFusionsTsv)) {
            LOGGER.warn(LINX_FUSION_TSV + " has to be an existing file: " + linxFusionsTsv);
        } else if (linxDisruptionsTsv == null || !exists(linxDisruptionsTsv)) {
            LOGGER.warn(LINX_DISRUPTION_TSV + " has to be an existing file: " + linxDisruptionsTsv);
        } else if (bachelorCsv != null && !exists(bachelorCsv)) {
            // Note: Bachelor CSV is optional (only exists in case pathogenic germline variants have been found).
            LOGGER.warn(BACHELOR_CSV + " has to be an existing file: " + bachelorCsv);
        } else if (chordPredictionFile == null || !exists(chordPredictionFile)) {
            LOGGER.warn(CHORD_PREDICTION_FILE + " has to be an existing file: " + chordPredictionFile);
        } else if (circosFile == null || !exists(circosFile)) {
            LOGGER.warn(CIRCOS_FILE + " has to be an existing file: " + circosFile);
        } else if (refSample == null) {
            LOGGER.warn(REF_SAMPLE + " has to be provided.");
        } else if (knowledgebaseDirectory == null || !exists(knowledgebaseDirectory) || !isDirectory(knowledgebaseDirectory)) {
            LOGGER.warn(KNOWLEDGEBASE_DIRECTORY + " has to be an existing directory: " + knowledgebaseDirectory);
        } else if (drupGenesCsv == null || !exists(drupGenesCsv)) {
            LOGGER.warn(DRUP_GENES_CSV + " has to be an existing file: " + drupGenesCsv);
        } else if (hotspotTsv == null || !exists(hotspotTsv)) {
            LOGGER.warn(HOTSPOT_TSV + " has to be an existing file: " + hotspotTsv);
        } else if (germlineGenesCsv == null || !exists(germlineGenesCsv)) {
            LOGGER.warn(GERMLINE_GENES_CSV + " has to be an existing file: " + germlineGenesCsv);
        } else if (sampleSummaryCsv == null || !exists(sampleSummaryCsv)) {
            LOGGER.warn(SAMPLE_SUMMARY_CSV + " has to be an existing file: " + sampleSummaryCsv);
        } else if (fastaFileLocation == null || !exists(fastaFileLocation)) {
            LOGGER.warn(FASTA_FILE_LOCATION + " has to be an existing file: " + fastaFileLocation);
        } else if (highConfidenceBed == null || !exists(highConfidenceBed)) {
            LOGGER.warn(HIGH_CONFIDENCE_BED + " has to be an existing file: " + highConfidenceBed);
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForQCFailReport(@NotNull final CommandLine cmd) {
        final QCFailReason qcFailReason = QCFailReason.fromIdentifier(cmd.getOptionValue(QC_FAIL_REASON));
        if (qcFailReason == QCFailReason.UNDEFINED) {
            LOGGER.warn(QC_FAIL_REASON + " has to be s, low_dna_yield, post_analysis_fail or shallow_seq.");
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForReportWriter(@NotNull final CommandLine cmd) {
        final String tumorSample = cmd.getOptionValue(TUMOR_SAMPLE);
        final String outputDirectory = cmd.getOptionValue(OUTPUT_DIRECTORY);

        if (tumorSample == null) {
            LOGGER.warn(TUMOR_SAMPLE + " has to be provided.");
        } else if (outputDirectory == null || !exists(outputDirectory) || !isDirectory(outputDirectory)) {
            LOGGER.warn(OUTPUT_DIRECTORY + " has to be an existing directory: " + outputDirectory);
        } else {
            return true;
        }
        return false;
    }

    private static boolean validInputForBaseReportData(@NotNull final CommandLine cmd) {
        final String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        final String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        final String hospitalDirectory = cmd.getOptionValue(HOSPITAL_DIRECTORY);

        final String signaturePath = cmd.getOptionValue(SIGNATURE);
        final String rvaLogoPath = cmd.getOptionValue(RVA_LOGO);
        final String companyLogoPath = cmd.getOptionValue(COMPANY_LOGO);

        if (tumorLocationCsv == null || !exists(tumorLocationCsv)) {
            LOGGER.warn(TUMOR_LOCATION_CSV + " has to be an existing file: " + tumorLocationCsv);
        } else if (limsDirectory == null || !exists(limsDirectory) || !isDirectory(limsDirectory)) {
            LOGGER.warn(LIMS_DIRECTORY + " has to be an existing directory: " + limsDirectory);
        } else if (hospitalDirectory == null || !exists(hospitalDirectory) || !isDirectory(hospitalDirectory)) {
            LOGGER.warn(HOSPITAL_DIRECTORY + " has to be an existing directory: " + hospitalDirectory);
        } else if (signaturePath == null || !exists(signaturePath)) {
            LOGGER.warn(SIGNATURE + " has to be an existing file: " + signaturePath);
        } else if (rvaLogoPath == null || !exists(rvaLogoPath)) {
            LOGGER.warn(RVA_LOGO + " has to be an existing file: " + rvaLogoPath);
        } else if (companyLogoPath == null || !exists(companyLogoPath)) {
            LOGGER.warn(COMPANY_LOGO + " has to be an existing file: " + companyLogoPath);
        } else {
            return true;
        }
        return false;
    }

    private static boolean exists(@NotNull final String path) {
        return Files.exists(new File(path).toPath());
    }

    private static boolean isDirectory(@NotNull final String path) {
        return Files.isDirectory(new File(path).toPath());
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(TUMOR_SAMPLE, true, "The sample for which a patient report will be generated.");
        options.addOption(OUTPUT_DIRECTORY, true, "Path to where the PDF reports have to be written to.");

        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the directory holding the LIMS data");
        options.addOption(HOSPITAL_DIRECTORY, true, "Path towards the directory containing hospital data.");

        options.addOption(RVA_LOGO, true, "Path towards a image file containing the RVA logo.");
        options.addOption(COMPANY_LOGO, true, "Path towards a image file containing the company logo.");
        options.addOption(SIGNATURE, true, "Path towards a image file containing the signature to be appended at the end of the report.");

        options.addOption(QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(QC_FAIL_REASON,
                true,
                "Either 'low_tumor_percentage', 'low_dna_yield', 'post_analysis_fail', 'shallow_seq' or 'insufficient_tissue_delivered'");

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(LINX_DISRUPTION_TSV, true, "Path towards the linx disruption TSV.");
        options.addOption(BACHELOR_CSV, true, "Path towards the germline germline CSV (optional).");
        options.addOption(CHORD_PREDICTION_FILE, true, "Path towards the CHORD prediction file.");

        options.addOption(REF_SAMPLE, true, "The reference sample for the sample for which we are generating a report.");
        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the directory holding knowledgebase output files.");
        options.addOption(FASTA_FILE_LOCATION, true, "Path towards the FASTA file containing the ref genome.");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path towards the high confidence BED file.");
        options.addOption(DRUP_GENES_CSV, true, "Path towards a CSV containing genes that could potentially indicate inclusion in DRUP.");
        options.addOption(HOTSPOT_TSV, true, "Path towards a TSV containing known hotspot variants.");
        options.addOption(GERMLINE_GENES_CSV, true, "Path towards a CSV containing germline genes which we want to report.");
        options.addOption(SAMPLE_SUMMARY_CSV, true, "Path towards a CSV containing the (clinical) summaries of the samples.");

        options.addOption(COMMENTS, true, "Additional comments to be added to the report (optional).");
        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Patient-Reporter", options);
        System.exit(1);
    }
}
