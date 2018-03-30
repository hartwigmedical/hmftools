package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusModel;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusReader;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.readers.PatientReader;
import com.hartwig.hmftools.patientdb.readers.RunsFolderReader;
import com.hartwig.hmftools.patientdb.validators.CurationValidator;
import com.hartwig.hmftools.patientdb.validators.PatientValidator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LoadClinicalData {
    private static final Logger LOGGER = LogManager.getLogger(LoadClinicalData.class);
    private static final String VERSION = LoadClinicalData.class.getPackage().getImplementationVersion();

    private static final String RUNS_DIR = "runs_dir";
    private static final String ECRF_FILE = "ecrf";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String LIMS_JSON = "lims_json";
    private static final String PRE_LIMS_ARRIVAL_DATES_CSV = "pre_lims_arrival_dates_csv";
    private static final String FORM_STATUS_CSV = "form_status_csv";
    private static final String DO_LOAD_RAW_ECRF = "do_load_raw_ecrf";
    private static final String CSV_OUT_DIR = "csv_out_dir";
    private static final String CANCER_TYPES_LINK = "cancer_types_symlink";
    private static final String PORTAL_DATA_LINK = "portal_data_symlink";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, XMLStreamException, SQLException {
        LOGGER.info("Running patient-db v{}", VERSION);
        final Options basicOptions = createBasicOptions();
        final Options clinicalOptions = createLimsOptions();
        final Options ecrfOptions = createEcrfOptions();
        final Options options = mergeOptions(basicOptions, clinicalOptions, ecrfOptions);

        final CommandLine cmd = createCommandLine(args, options);
        final String runsFolderPath = cmd.getOptionValue(RUNS_DIR);
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);

        final boolean loadRawEcrf = cmd.hasOption(DO_LOAD_RAW_ECRF);

        if (Utils.anyNull(runsFolderPath, userName, password, databaseUrl)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db", options);
        } else {
            final File runDirectory = new File(runsFolderPath);
            if (runDirectory.isDirectory()) {
                LOGGER.info("Running clinical data import.");
                final List<RunContext> runContexts = RunsFolderReader.getRunContexts(runDirectory);
                final String jdbcUrl = "jdbc:" + databaseUrl;

                final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);
                if (loadRawEcrf) {
                    writeRawEcrf(ecrfOptions, cmd, runContexts, dbWriter);
                }

                writeClinicalData(clinicalOptions, cmd, runContexts, dbWriter);
            } else {
                if (!runDirectory.exists()) {
                    LOGGER.warn("dir " + runDirectory + " does not exist.");
                }
                final HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("patient-db", basicOptions);
            }
        }
    }

    private static void writeClinicalData(@NotNull final Options clinicalOptions, @NotNull final CommandLine cmd,
            @NotNull final List<RunContext> runContexts, @NotNull final DatabaseAccess dbAccess) throws IOException, XMLStreamException {
        final String ecrfFilePath = cmd.getOptionValue(ECRF_FILE);
        final String limsJson = cmd.getOptionValue(LIMS_JSON);
        final String preLIMSArrivalDatesCsv = cmd.getOptionValue(PRE_LIMS_ARRIVAL_DATES_CSV);
        final String formStatusCsv = cmd.getOptionValue(FORM_STATUS_CSV);
        final String csvOutputDir = cmd.getOptionValue(CSV_OUT_DIR);
        final Optional<String> cancerTypesLink = Optional.ofNullable(cmd.getOptionValue(CANCER_TYPES_LINK));
        final Optional<String> portalDataLink = Optional.ofNullable(cmd.getOptionValue(PORTAL_DATA_LINK));

        if (Utils.anyNull(ecrfFilePath, limsJson, preLIMSArrivalDatesCsv, formStatusCsv, csvOutputDir)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db", clinicalOptions);
        } else {
            LOGGER.info("Loading ecrf model...");
            dbAccess.clearClinicalTables();

            FormStatusModel formStatusModel = FormStatusReader.buildModelFromCsv(formStatusCsv);
            CpctEcrfModel model = CpctEcrfModel.loadFromXML(ecrfFilePath, formStatusModel);
            Lims lims = LimsFactory.fromLimsJsonWithPreLIMSArrivalDates(limsJson, preLIMSArrivalDatesCsv);
            TreatmentCurator treatmentCurator = TreatmentCurator.fromProductionResource();
            TumorLocationCurator tumorLocationCurator = TumorLocationCurator.fromProductionResource();
            final PatientReader patientReader = new PatientReader(model, treatmentCurator, tumorLocationCurator, lims);

            final Set<String> sequencedCpctPatientIds = Utils.sequencedPatientIds(runContexts)
                    .stream()
                    .filter(patientId -> patientId.startsWith("CPCT"))
                    .collect(Collectors.toSet());

            final Map<String, Patient> readPatients = readEcrfPatients(patientReader, model.patients(), runContexts);
            LOGGER.info("Read {} patients from ecrf", readPatients.size());
            DumpClinicalData.writeClinicalDumps(csvOutputDir, readPatients.values(), cancerTypesLink, portalDataLink);

            LOGGER.info("Writing CPCT clinical data for {} sequenced patients.", sequencedCpctPatientIds.size());
            for (final String patientId : sequencedCpctPatientIds) {
                final Patient patient = readPatients.get(patientId);
                if (patient == null) {
                    LOGGER.error("Could not find patient with id: " + patientId + " in ecrf file.");
                } else {
                    dbAccess.writeClinicalData(patient);
                    final List<ValidationFinding> findings = PatientValidator.validatePatient(patient);

                    dbAccess.writeValidationFindings(findings);
                    dbAccess.writeValidationFindings(patient.matchFindings());
                }
            }
            dbAccess.writeValidationFindings(CurationValidator.validateTreatmentCurator(treatmentCurator));
            dbAccess.writeValidationFindings(CurationValidator.validateTumorLocationCurator(tumorLocationCurator));

            LOGGER.info("Done!");
        }
    }

    @NotNull
    private static Map<String, Patient> readEcrfPatients(@NotNull final PatientReader reader, @NotNull final Iterable<EcrfPatient> patients,
            @NotNull final List<RunContext> runContexts) throws IOException {
        final Map<String, Patient> patientMap = Maps.newHashMap();
        for (final EcrfPatient ecrfPatient : patients) {
            final List<String> tumorSamplesForPatient = getTumorSamplesForPatient(ecrfPatient.patientId(), runContexts);
            final Patient patient = reader.read(ecrfPatient, tumorSamplesForPatient);
            patientMap.put(patient.patientData().cpctId(), patient);
        }
        return patientMap;
    }

    private static void writeRawEcrf(@NotNull final Options ecrfOptions, @NotNull final CommandLine cmd,
            @NotNull final List<RunContext> runContexts, @NotNull final DatabaseAccess dbWriter) throws IOException, XMLStreamException {
        final String ecrfFilePath = cmd.getOptionValue(ECRF_FILE);
        final String formStatusPath = cmd.getOptionValue(FORM_STATUS_CSV);
        if (Utils.anyNull(ecrfFilePath, formStatusPath)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db -" + DO_LOAD_RAW_ECRF, ecrfOptions);
        } else {
            dbWriter.clearCpctEcrf();
            LOGGER.info("Loading ecrf model...");
            final FormStatusModel formStatusModel = FormStatusReader.buildModelFromCsv(formStatusPath);
            final CpctEcrfModel model = CpctEcrfModel.loadFromXML(ecrfFilePath, formStatusModel);
            final Set<String> sequencedPatients = Utils.sequencedPatientIds(runContexts);
            LOGGER.info("Writing raw ecrf data for " + sequencedPatients.size() + " patients.");
            dbWriter.writeCpctEcrf(model, sequencedPatients);
            LOGGER.info("Done writing raw ecrf data for " + sequencedPatients.size() + " patients!");
        }
    }

    @NotNull
    private static List<String> getTumorSamplesForPatient(@NotNull final String patientId, @NotNull final List<RunContext> runContexts) {
        final List<String> sampleIdsForPatient = Lists.newArrayList();
        runContexts.forEach(runContext -> {
            final String sampleId = runContext.tumorSample();
            if (sampleId.startsWith(patientId) && !sampleIdsForPatient.contains(sampleId)) {
                sampleIdsForPatient.add(sampleId);
            }
        });
        return sampleIdsForPatient;
    }

    @NotNull
    private static Options mergeOptions(@NotNull final Options... optionsArray) {
        final Options options = new Options();
        final List<Options> optionsList = Lists.newArrayList(optionsArray);
        optionsList.forEach(opt -> opt.getOptions().forEach(options::addOption));
        return options;
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIR, true, "Path towards the folder containing patient runs.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(DO_LOAD_RAW_ECRF, false, "Also write raw ecrf data to database?");
        return options;
    }

    @NotNull
    private static Options createLimsOptions() {
        final Options options = new Options();
        options.addOption(LIMS_JSON, true, "Path towards the LIMS json file.");
        options.addOption(PRE_LIMS_ARRIVAL_DATES_CSV, true, "Path towards the pre-HMF arrival date csv.");
        return options;
    }

    @NotNull
    private static Options createEcrfOptions() {
        final Options options = new Options();
        options.addOption(ECRF_FILE, true, "Path towards the cpct ecrf file.");
        options.addOption(FORM_STATUS_CSV, true, "Path towards the form status csv file.");
        options.addOption(CSV_OUT_DIR, true, "Path towards the output directory for csv data dumps.");
        options.addOption(CANCER_TYPES_LINK, true, "Name of cancer type csv symlink.");
        options.addOption(PORTAL_DATA_LINK, true, "Name of portal data csv symlink.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
