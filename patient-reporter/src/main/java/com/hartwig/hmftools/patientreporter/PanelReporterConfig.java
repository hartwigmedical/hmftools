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
public interface PanelReporterConfig {

    Logger LOGGER = LogManager.getLogger(PanelReporterConfig.class);

    // General params needed for every report
    String TUMOR_SAMPLE_ID = "tumor_sample_id";
    String TUMOR_SAMPLE_BARCODE = "tumor_sample_barcode";
    String OUTPUT_DIRECTORY_REPORT = "output_dir_report";

    String PRIMARY_TUMOR_TSV = "primary_tumor_tsv";
    String LIMS_DIRECTORY = "lims_dir";

    String COMPANY_LOGO = "company_logo";
    String SIGNATURE = "signature";

    // General params needed for every report but for QC fail it can be optional in some cases
    String REF_SAMPLE_ID = "ref_sample_id";
    String REF_SAMPLE_BARCODE = "ref_sample_barcode";

    // Params specific for Panel reports
    String PANEL = "panel";
    String PANEL_QC_FAIL = "panel_qc_fail";
    String PANEL_QC_FAIL_REASON = "panel_qc_fail_reason";
    String PANEL_VCF_NAME = "panel_vcf_name";

    // Some additional optional params and flags
    String COMMENTS = "comments";
    String CORRECTED_REPORT = "corrected_report";
    String CORRECTED_REPORT_EXTERN = "corrected_report_extern";
    String LOG_DEBUG = "log_debug";

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

        options.addOption(PRIMARY_TUMOR_TSV, true, "Path towards the (curated) primary tumor TSV.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the directory holding the LIMS data");

        options.addOption(COMPANY_LOGO, true, "Path towards an image file containing the company logo.");
        options.addOption(SIGNATURE, true, "Path towards an image file containing the signature to be appended at the end of the report.");

        options.addOption(REF_SAMPLE_ID, true, "The reference sample ID for the tumor sample for which a report is generated.");
        options.addOption(REF_SAMPLE_BARCODE, true, "The reference sample barcode for the tumor sample for which a report is generated.");

        options.addOption(PANEL, false, "If set, generates a panel report.");
        options.addOption(PANEL_QC_FAIL, false, "If set, generates a qc-fail report.");
        options.addOption(PANEL_QC_FAIL_REASON, true, "One of: " + Strings.join(Lists.newArrayList(QCFailReason.validIdentifiers()), ','));
        options.addOption(PANEL_VCF_NAME, true, "The name of the VCF file of the panel results.");

        options.addOption(COMMENTS, true, "Additional comments to be added to the report (optional).");
        options.addOption(CORRECTED_REPORT, false, "If provided, generate a corrected report with corrected name");
        options.addOption(CORRECTED_REPORT_EXTERN, false, "If provided, generate a corrected report with intern/extern correction");

        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");

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
    String primaryTumorTsv();

    @NotNull
    String limsDir();

    @NotNull
    String companyLogo();

    @NotNull
    String signature();

    boolean panel();

    boolean panelQcFail();

    @NotNull
    String panelVCFname();

    @Nullable
    PanelFailReason panelQcFailReason();

    @Nullable
    String comments();

    boolean isCorrectedReport();

    boolean isCorrectedReportExtern();

    boolean requirePipelineVersionFile();

    @Nullable
    String pipelineVersionFile();

    @NotNull
    String expectedPipelineVersion();

    boolean overridePipelineVersion();

    @NotNull
    RefGenomeVersion refGenomeVersion();

    @NotNull
    static PanelReporterConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
            LOGGER.debug("Switched root level logging to DEBUG");
        }

        boolean requirePipelineVersion = cmd.hasOption(REQUIRE_PIPELINE_VERSION_FILE);

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

            panelVCFFile = nonOptionalValue(cmd, PANEL_VCF_NAME);
        }

        return ImmutablePanelReporterConfig.builder()
                .refSampleId(cmd.hasOption(REF_SAMPLE_ID) ? nonOptionalValue(cmd, REF_SAMPLE_ID) : null)
                .refSampleBarcode(cmd.hasOption(REF_SAMPLE_BARCODE) ? nonOptionalValue(cmd, REF_SAMPLE_BARCODE) : null)
                .tumorSampleId(nonOptionalValue(cmd, TUMOR_SAMPLE_ID))
                .tumorSampleBarcode(nonOptionalValue(cmd, TUMOR_SAMPLE_BARCODE))
                .outputDirReport(nonOptionalDir(cmd, OUTPUT_DIRECTORY_REPORT))
                .primaryTumorTsv(nonOptionalFile(cmd, PRIMARY_TUMOR_TSV))
                .limsDir(nonOptionalDir(cmd, LIMS_DIRECTORY))
                .companyLogo(nonOptionalFile(cmd, COMPANY_LOGO))
                .signature(nonOptionalFile(cmd, SIGNATURE))
                .panel(isPanel)
                .panelQcFail(isPanelQCFail)
                .panelQcFailReason(panelQcFailReason)
                .panelVCFname(panelVCFFile)
                .comments(cmd.getOptionValue(COMMENTS))
                .isCorrectedReport(cmd.hasOption(CORRECTED_REPORT))
                .isCorrectedReportExtern(cmd.hasOption(CORRECTED_REPORT_EXTERN))
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
