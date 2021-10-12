package com.hartwig.hmftools.patientdb.clinical;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.File;
import java.nio.file.Files;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ClinicalAlgoConfig {

    String RUNS_DIRECTORY = "runs_dir";
    String RUNS_JSON = "runs_json";
    String PIPELINE_VERSION_FILE = "pipeline_version_file";

    String CPCT_ECRF_FILE = "cpct_ecrf";
    String CPCT_FORM_STATUS_CSV = "cpct_form_status_csv";
    String DRUP_ECRF_FILE = "drup_ecrf";

    String CONSENT_CONFIG_TSV = "consent_config_tsv";

    String DO_LOAD_CLINICAL_DATA = "do_load_clinical_data";
    String DO_LOAD_RAW_ECRF = "do_load_raw_ecrf";

    String DO_PROCESS_WIDE_CLINICAL_DATA = "do_process_wide_clinical_data";
    String WIDE_PRE_AVL_TREATMENT_CSV = "wide_pre_avl_treatment_csv";
    String WIDE_BIOPSY_CSV = "wide_biopsy_csv";
    String WIDE_AVL_TREATMENT_CSV = "wide_avl_treatment_csv";
    String WIDE_RESPONSE_CSV = "wide_response_csv";
    String WIDE_FIVE_DAYS_CSV = "wide_five_days_csv";

    String LIMS_DIRECTORY = "lims_dir";
    String REPORTING_DB_TSV = "reporting_db_tsv";

    String CURATED_PRIMARY_TUMOR_TSV = "curated_primary_tumor_tsv";
    String PATIENT_TUMOR_CURATION_STATUS_TSV = "patient_tumor_curation_status_tsv";

    String DOID_JSON = "doid_json";
    String TUMOR_LOCATION_MAPPING_TSV = "tumor_location_mapping_tsv";
    String TUMOR_LOCATION_OVERRIDES_TSV = "tumor_location_overrides_tsv";
    String BIOPSY_MAPPING_TSV = "biopsy_mapping_tsv";
    String TREATMENT_MAPPING_TSV = "treatment_mapping_tsv";

    @NotNull
    static Options createOptions() {
        Options options = new Options();
        options.addOption(RUNS_DIRECTORY,
                true,
                "Path towards the folder containing patient runs that are considered part of HMF database.");
        options.addOption(RUNS_JSON, true, "Path towards a JSON file containing patient runs that are considered part of HMF database.");
        options.addOption(PIPELINE_VERSION_FILE, true, "Path towards the pipeline version");

        options.addOption(CPCT_ECRF_FILE, true, "Path towards the CPCT ecrf file.");
        options.addOption(CPCT_FORM_STATUS_CSV, true, "Path towards the CPCT form status csv file.");
        options.addOption(DRUP_ECRF_FILE, true, "Path towards the DRUP ecrf file.");

        options.addOption(CONSENT_CONFIG_TSV, true, "Path towards the informed consent config TSV file.");

        options.addOption(DO_LOAD_CLINICAL_DATA, false, "If set, clinical data will be loaded into the database.");
        options.addOption(DO_LOAD_RAW_ECRF, false, "If set, writes raw ecrf data to database.");

        options.addOption(DO_PROCESS_WIDE_CLINICAL_DATA,
                false,
                "if set, creates clinical timeline for WIDE patients and persists to database.");
        options.addOption(WIDE_AVL_TREATMENT_CSV, true, "Path towards the WIDE avl treatment csv");
        options.addOption(WIDE_PRE_AVL_TREATMENT_CSV, true, "Path towards the WIDE pre avl treatment csv.");
        options.addOption(WIDE_BIOPSY_CSV, true, "Path towards the WIDE biopsy csv.");
        options.addOption(WIDE_RESPONSE_CSV, true, "Path towards the WIDE response csv.");
        options.addOption(WIDE_FIVE_DAYS_CSV, true, "Path towards the WIDE five days csv.");

        options.addOption(LIMS_DIRECTORY, true, "Path towards the LIMS directory.");
        options.addOption(REPORTING_DB_TSV, true, "Path towards the reporting db tsv file.");

        options.addOption(CURATED_PRIMARY_TUMOR_TSV, true, "Path towards to the curated primary tumor TSV.");
        options.addOption(PATIENT_TUMOR_CURATION_STATUS_TSV, true, "Path where patient tumor curation status will be written to");

        options.addOption(DOID_JSON, true, "Path towards to the json file of the doid ID of primary tumors.");
        options.addOption(TUMOR_LOCATION_MAPPING_TSV, true, "Path towards to the tumor location mapping TSV.");
        options.addOption(TUMOR_LOCATION_OVERRIDES_TSV, true, "Path towards to the tumor location overrides TSV.");
        options.addOption(BIOPSY_MAPPING_TSV, true, "Path towards to the biopsy mapping TSV.");
        options.addOption(TREATMENT_MAPPING_TSV, true, "Path towards to the treatment mapping TSV.");

        addDatabaseCmdLineArgs(options);
        return options;
    }

    @Nullable
    String runsDirectory();

    @Nullable
    String runsJson();

    @Nullable
    String pipelineVersionFile();

    @NotNull
    String cpctEcrfFile();

    @NotNull
    String cpctFormStatusCsv();

    @NotNull
    String drupEcrfFile();

    @NotNull
    String consentConfigTsv();

    boolean doLoadClinicalData();

    boolean doLoadRawEcrf();

    boolean doProcessWideClinicalData();

    @Nullable
    String widePreAvlTreatmentCsv();

    @Nullable
    String wideBiopsyCsv();

    @Nullable
    String wideAvlTreatmentCsv();

    @Nullable
    String wideResponseCsv();

    @Nullable
    String wideFiveDaysCsv();

    @NotNull
    String limsDirectory();

    @NotNull
    String reportingDbTsv();

    @NotNull
    String curatedPrimaryTumorTsv();

    @NotNull
    String patientTumorCurationStatusTsv();

    @NotNull
    String doidJson();

    @NotNull
    String tumorLocationMappingTsv();

    @NotNull
    String tumorLocationOverridesTsv();

    @NotNull
    String biopsyMappingTsv();

    @NotNull
    String treatmentMappingTsv();

    @NotNull
    static ClinicalAlgoConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        String runsJson = optionalFile(cmd, RUNS_JSON);
        String runsDirectory = optionalDir(cmd, RUNS_DIRECTORY);

        if ((runsJson == null && runsDirectory == null) || (runsJson != null && runsDirectory != null)) {
            throw new IllegalStateException("Either a runs directory or runs json must be specified, and not both");
        }

        boolean doProcessWideClinicalData = cmd.hasOption(DO_PROCESS_WIDE_CLINICAL_DATA);
        ImmutableClinicalAlgoConfig.Builder builder = ImmutableClinicalAlgoConfig.builder()
                .runsDirectory(runsDirectory)
                .runsJson(runsJson)
                .pipelineVersionFile(nonOptionalValue(cmd, PIPELINE_VERSION_FILE))
                .cpctEcrfFile(nonOptionalFile(cmd, CPCT_ECRF_FILE))
                .cpctFormStatusCsv(nonOptionalFile(cmd, CPCT_FORM_STATUS_CSV))
                .drupEcrfFile(nonOptionalFile(cmd, DRUP_ECRF_FILE))
                .consentConfigTsv(nonOptionalFile(cmd, CONSENT_CONFIG_TSV))
                .doLoadClinicalData(cmd.hasOption(DO_LOAD_CLINICAL_DATA))
                .doLoadRawEcrf(cmd.hasOption(DO_LOAD_RAW_ECRF))
                .doProcessWideClinicalData(doProcessWideClinicalData)
                .limsDirectory(nonOptionalDir(cmd, LIMS_DIRECTORY))
                .reportingDbTsv(nonOptionalFile(cmd, REPORTING_DB_TSV))
                .curatedPrimaryTumorTsv(nonOptionalValue(cmd, CURATED_PRIMARY_TUMOR_TSV))
                .patientTumorCurationStatusTsv(nonOptionalValue(cmd, PATIENT_TUMOR_CURATION_STATUS_TSV))
                .doidJson(nonOptionalFile(cmd, DOID_JSON))
                .tumorLocationMappingTsv(nonOptionalFile(cmd, TUMOR_LOCATION_MAPPING_TSV))
                .tumorLocationOverridesTsv(nonOptionalFile(cmd, TUMOR_LOCATION_OVERRIDES_TSV))
                .biopsyMappingTsv(nonOptionalFile(cmd, BIOPSY_MAPPING_TSV))
                .treatmentMappingTsv(nonOptionalFile(cmd, TREATMENT_MAPPING_TSV));

        if (doProcessWideClinicalData) {
            builder = builder.widePreAvlTreatmentCsv(nonOptionalFile(cmd, WIDE_PRE_AVL_TREATMENT_CSV))
                    .wideBiopsyCsv(nonOptionalFile(cmd, WIDE_BIOPSY_CSV))
                    .wideAvlTreatmentCsv(nonOptionalFile(cmd, WIDE_AVL_TREATMENT_CSV))
                    .wideResponseCsv(nonOptionalFile(cmd, WIDE_RESPONSE_CSV))
                    .wideFiveDaysCsv(nonOptionalFile(cmd, WIDE_FIVE_DAYS_CSV));
        }

        return builder.build();
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

    @Nullable
    static String optionalDir(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = cmd.getOptionValue(param);

        if (value == null) {
            return null;
        }

        if (!pathExists(value) || !pathIsDirectory(value)) {
            throw new ParseException("Parameter '" + param + "', if provided, must be an existing directory: " + value);
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

    @Nullable
    static String optionalFile(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = cmd.getOptionValue(param);

        if (value == null) {
            return null;
        }

        if (!pathExists(value)) {
            throw new ParseException("Parameter '" + param + "', if provided, must be an existing file: " + value);
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
