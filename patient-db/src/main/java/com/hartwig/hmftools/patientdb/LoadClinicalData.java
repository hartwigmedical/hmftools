package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusModel;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusReader;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.lims.LimsSampleType;
import com.hartwig.hmftools.patientdb.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.readers.ColoPatientReader;
import com.hartwig.hmftools.patientdb.readers.EcrfPatientReader;
import com.hartwig.hmftools.patientdb.readers.LimsPatientReader;
import com.hartwig.hmftools.patientdb.readers.LimsSampleReader;
import com.hartwig.hmftools.patientdb.readers.RunsFolderReader;
import com.hartwig.hmftools.patientdb.readers.cpct.CpctPatientReader;
import com.hartwig.hmftools.patientdb.readers.cpct.CpctUtil;
import com.hartwig.hmftools.patientdb.readers.drup.DrupPatientReader;
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

    private static final String RUNS_DIR_DATABASE = "runs_dir_db";
    private static final String RUNS_DIR_NON_DATABASE = "runs_dir_non_db";
    private static final String CPCT_ECRF_FILE = "cpct_ecrf";
    private static final String CPCT_FORM_STATUS_CSV = "cpct_form_status_csv";
    private static final String DRUP_ECRF_FILE = "drup_ecrf";
    private static final String DO_LOAD_RAW_ECRF = "do_load_raw_ecrf";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String LIMS_DIRECTORY = "lims";

    private static final String CSV_OUT_DIR = "csv_out_dir";
    private static final String TUMOR_LOCATION_SYMLINK = "tumor_location_symlink";
    private static final String PORTAL_DATA_LINK = "portal_data_symlink";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, XMLStreamException, SQLException {
        LOGGER.info("Running patient-db v{}", VERSION);
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (checkInputs(cmd)) {
            LOGGER.info("Running clinical data import.");

            final DatabaseAccess dbWriter = createDbWriter(cmd);

            final Lims lims = buildLims(cmd);

            final Map<String, List<SampleData>> samplesPerPatientAll = loadSamplesPerPatientAll(lims);
            final Map<String, List<SampleData>> samplesPerPatientFilteredToOnlySequenced = filterPatients(cmd, samplesPerPatientAll);
            final EcrfModels ecrfModels = loadEcrfModels(cmd);

            if (cmd.hasOption(DO_LOAD_RAW_ECRF)) {
                writeRawEcrf(dbWriter, samplesPerPatientFilteredToOnlySequenced.keySet(), ecrfModels);
            }

            writeClinicalData(dbWriter,
                    samplesPerPatientFilteredToOnlySequenced,
                    samplesPerPatientAll,
                    ecrfModels,
                    cmd.getOptionValue(CSV_OUT_DIR),
                    Optional.ofNullable(cmd.getOptionValue(TUMOR_LOCATION_SYMLINK)),
                    Optional.ofNullable(cmd.getOptionValue(PORTAL_DATA_LINK)));
        } else {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db", options);
        }
    }

    @NotNull
    private static Map<String, List<SampleData>> filterPatients(@NotNull CommandLine cmd,
            @NotNull Map<String, List<SampleData>> samplesPerPatientAll) throws IOException {
        List<RunContext> contextsRunExtracted = loadRunContexts(cmd);

        Map<String, List<SampleData>> samplesPerPatientSequenced = extractSamplesFromRunContexts(contextsRunExtracted, samplesPerPatientAll);
        LOGGER.info(String.format("Using sequenced samples for %s patients from LIMS", samplesPerPatientSequenced.keySet().size()));

        return samplesPerPatientSequenced;
    }

    @NotNull
    private static Map<String, List<SampleData>> extractSamplesFromRunContexts(@NotNull List<RunContext> runContexts,
            Map<String, List<SampleData>> samplesPerPatientAll) {

        final Set<String> sequencedPatientIdentifiers = Utils.sequencedPatientIdentifiers(runContexts);

        Map<String, List<SampleData>> samplesPerPatientSequenced = Maps.newHashMap();
        for (String patientIdentifier : sequencedPatientIdentifiers) {
            boolean filteredSequencedSamples = samplesPerPatientAll.keySet().contains(patientIdentifier);
            if (filteredSequencedSamples) {
                Set<String> sampleIds = extractTumorSampleIdsForPatient(patientIdentifier, runContexts);
                samplesPerPatientSequenced.put(patientIdentifier, samplesPerPatientAll.get(patientIdentifier));
                LOGGER.info(samplesPerPatientSequenced);
            }
        }
        return samplesPerPatientSequenced;
    }

    @NotNull
    private static Lims buildLims(@NotNull CommandLine cmd) throws IOException {
        final String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        LOGGER.info(String.format("Loading samples from LIMS on %s.", limsDirectory));
        return LimsFactory.fromLimsDirectory(limsDirectory);
    }

    private static void writeClinicalData(@NotNull final DatabaseAccess dbAccess,
            @NotNull Map<String, List<SampleData>> samplesPerPatientSequenced, Map<String, List<SampleData>> samplesPerPatientAll,
            @NotNull EcrfModels ecrfModels, @NotNull String csvOutputDir, @NotNull Optional<String> tumorLocationSymlink,
            @NotNull Optional<String> portalDataLink) throws IOException {
        TumorLocationCurator tumorLocationCurator = TumorLocationCurator.fromProductionResource();
        BiopsySiteCurator biopsySiteCurator = BiopsySiteCurator.fromProductionResource();
        TreatmentCurator treatmentCurator = TreatmentCurator.fromProductionResource();

        Map<String, Patient> patients =
                loadAndInterpretPatients(samplesPerPatientSequenced, ecrfModels, tumorLocationCurator, treatmentCurator, biopsySiteCurator);

        Map<String, Patient> patientsAll =
                loadAndInterpretPatients(samplesPerPatientAll, ecrfModels, tumorLocationCurator, treatmentCurator, biopsySiteCurator);

        DumpClinicalData.writeClinicalDumps(csvOutputDir, patientsAll.values(), tumorLocationSymlink, portalDataLink);

        LOGGER.info("Clearing interpreted clinical tables in database.");
        dbAccess.clearClinicalTables();

        Set<String> sequencedPatientIdentifiers = samplesPerPatientSequenced.keySet();
        int missingPatients = 0;
        int missingSamples = 0;
        LOGGER.info(String.format("Writing clinical data for %s sequenced patients.", sequencedPatientIdentifiers.size()));
        for (final String patientIdentifier : sequencedPatientIdentifiers) {
            Patient patient = patients.get(patientIdentifier);
            if (patient == null) {
                missingPatients++;
                missingSamples += samplesPerPatientSequenced.get(patientIdentifier).size();
                dbAccess.writeSampleClinicalData(patientIdentifier, samplesPerPatientSequenced.get(patientIdentifier));
            } else {
                dbAccess.writeFullClinicalData(patient);
                List<ValidationFinding> findings = PatientValidator.validatePatient(patient);

                dbAccess.writeValidationFindings(findings);
                dbAccess.writeValidationFindings(patient.matchFindings());
            }
        }
        if (missingPatients > 0) {
            LOGGER.warn(String.format("Could not load %s patients (%s samples)!", missingPatients, missingSamples));
        }
        dbAccess.writeValidationFindings(CurationValidator.validateTreatmentCurator(treatmentCurator));
        dbAccess.writeValidationFindings(CurationValidator.validateTumorLocationCurator(tumorLocationCurator));

        LOGGER.info("Finished!");
    }

    @NotNull
    private static Map<String, Patient> loadAndInterpretPatients(@NotNull Map<String, List<SampleData>> samplesPerPatient,
            @NotNull EcrfModels ecrfModels, @NotNull TumorLocationCurator tumorLocationCurator, @NotNull TreatmentCurator treatmentCurator,
            @NotNull BiopsySiteCurator biopsySiteCurator) {
        final EcrfModel cpctEcrfModel = ecrfModels.cpctModel();
        LOGGER.info(String.format("Interpreting and curating data for %s CPCT patients.", cpctEcrfModel.patientCount()));
        EcrfPatientReader cpctPatientReader = new CpctPatientReader(tumorLocationCurator,
                CpctUtil.extractHospitalMap(cpctEcrfModel),
                biopsySiteCurator,
                treatmentCurator);

        Map<String, Patient> cpctPatients = readEcrfPatients(cpctPatientReader, cpctEcrfModel.patients(), samplesPerPatient);
        LOGGER.info(String.format("Finished curation of %s CPCT patients.", cpctPatients.size()));

        final EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info(String.format("Interpreting and curating data for %s DRUP patients.", drupEcrfModel.patientCount()));
        EcrfPatientReader drupPatientReader = new DrupPatientReader(tumorLocationCurator, biopsySiteCurator);

        Map<String, Patient> drupPatients = readEcrfPatients(drupPatientReader, drupEcrfModel.patients(), samplesPerPatient);
        LOGGER.info(String.format("Finished curation of %s DRUP patients.", drupPatients.size()));

        LOGGER.info("Interpreting and curating data based off LIMS");
        Map<String, Patient> patientsFromLims = readLimsPatients(samplesPerPatient, tumorLocationCurator);
        LOGGER.info(String.format("Finished curation of %s patients based off LIMS", patientsFromLims.keySet().size()));

        Map<String, Patient> mergedPatients = Maps.newHashMap();
        mergedPatients.putAll(cpctPatients);
        mergedPatients.putAll(drupPatients);
        mergedPatients.putAll(patientsFromLims);
        mergedPatients.putAll(readColoPatients());
        return mergedPatients;
    }

    @NotNull
    private static Map<String, Patient> readEcrfPatients(@NotNull final EcrfPatientReader reader,
            @NotNull final Iterable<EcrfPatient> patients, @NotNull final Map<String, List<SampleData>> samplesPerPatientSequenced) {
        final Map<String, Patient> patientMap = Maps.newHashMap();
        for (final EcrfPatient ecrfPatient : patients) {
            List<SampleData> samples = samplesPerPatientSequenced.get(ecrfPatient.patientId());
            Patient patient = reader.read(ecrfPatient, samples != null ? samples : Lists.newArrayList());
            patientMap.put(patient.patientIdentifier(), patient);
        }
        return patientMap;
    }

    @NotNull
    private static Map<String, Patient> readLimsPatients(@NotNull final Map<String, List<SampleData>> samplesPerPatientSequenced,
            @NotNull TumorLocationCurator tumorLocationCurator) {
        final Map<String, Patient> patientMap = Maps.newHashMap();
        final LimsPatientReader limsPatientReader = new LimsPatientReader(tumorLocationCurator);

        for (Map.Entry<String, List<SampleData>> entry : samplesPerPatientSequenced.entrySet()) {
            List<SampleData> samples = entry.getValue();
            if (!samples.isEmpty()) {
                LimsSampleType sampleType = LimsSampleType.fromSampleId(samples.get(0).sampleId());

                if (sampleType == LimsSampleType.CORE || sampleType == LimsSampleType.WIDE) {
                    String patientId = entry.getKey();
                    Patient limsPatient = limsPatientReader.read(patientId, samples);
                    patientMap.put(patientId, limsPatient);
                }
            }
        }

        return patientMap;
    }

    @NotNull
    private static Map<String, Patient> readColoPatients() {
        final Map<String, Patient> patientMap = Maps.newHashMap();
        final ColoPatientReader coloPatientReader = new ColoPatientReader();
        LOGGER.info("Creating patient representation for COLO829");
        Patient colo829Patient = coloPatientReader.read("COLO829T");

        patientMap.put(colo829Patient.patientIdentifier(), colo829Patient);
        return patientMap;
    }

    private static void writeRawEcrf(@NotNull DatabaseAccess dbWriter, @NotNull Set<String> sequencedPatients,
            @NotNull EcrfModels ecrfModels) {
        final EcrfModel cpctEcrfModel = ecrfModels.cpctModel();
        LOGGER.info(String.format("Writing raw cpct ecrf data for %s patients", cpctEcrfModel.patientCount()));
        dbWriter.clearCpctEcrf();
        dbWriter.writeCpctEcrf(cpctEcrfModel, sequencedPatients);
        LOGGER.info(String.format("Finished writing raw cpct ecrf data for %s patients.", cpctEcrfModel.patientCount()));

        final EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info(String.format("Writing raw drup ecrf data for %s patients", drupEcrfModel.patientCount()));
        dbWriter.clearDrupEcrf();
        dbWriter.writeDrupEcrf(drupEcrfModel, sequencedPatients);
        LOGGER.info(String.format("Finished writing raw drup ecrf data for %s patients.", drupEcrfModel.patientCount()));
    }

    @NotNull
    private static List<RunContext> loadRunContexts(@NotNull CommandLine cmd) throws IOException {
        final String runsFolderPathDb = cmd.getOptionValue(RUNS_DIR_DATABASE);
        final List<RunContext> runContextsDb = RunsFolderReader.getRunContexts(new File(runsFolderPathDb));
        LOGGER.info(String.format("Loading run contexts from %s (%s sets)", runsFolderPathDb, runContextsDb.size()));

        final String runsFolderPathNonDb = cmd.getOptionValue(RUNS_DIR_NON_DATABASE);
        final List<RunContext> runContextsNonDb = RunsFolderReader.getRunContexts(new File(runsFolderPathNonDb));
        LOGGER.info(String.format("Loading run contexts from %s (%s sets)", runsFolderPathNonDb, runContextsNonDb.size()));

        List<RunContext> runContextsAll = Lists.newArrayList();
        runContextsAll.addAll(runContextsDb);
        runContextsAll.addAll(runContextsNonDb);

        LOGGER.info(String.format("Finished loading %s run contexts.", runContextsAll.size()));
        return runContextsAll;
    }

    @NotNull
    private static Map<String, List<SampleData>> loadSamplesPerPatientAll(@NotNull Lims lims) {
        Map<String, List<SampleData>> samplesPerPatientAll = extractSamplesFromLims(lims);
        LOGGER.info(String.format("Loaded all samples for %s patients from LIMS", samplesPerPatientAll.keySet().size()));

        return samplesPerPatientAll;
    }

    @NotNull
    private static Map<String, List<SampleData>> extractSamplesFromLims(@NotNull Lims lims) {
        LimsSampleReader sampleReader = new LimsSampleReader(lims);
        Set<String> allSampleIDs = lims.allSampleIds();

        Map<String, List<SampleData>> samplesPerPatientAll = Maps.newHashMap();
        for (String sampleId : allSampleIDs) {
            LimsSampleType sampleType = LimsSampleType.fromSampleId(sampleId);
            if (sampleType.equals(LimsSampleType.CORE) || sampleType.equals(LimsSampleType.WIDE) || sampleType.equals(LimsSampleType.DRUP)
                    || sampleType.equals(LimsSampleType.CPCT)) {
                if (sampleId.substring(12).contains("T")) {
                    String patientIdentifier = sampleId.substring(0, 12);
                    samplesPerPatientAll.put(patientIdentifier, sampleReader.read(Sets.newHashSet(sampleId)));
                }
            }
        }
        return samplesPerPatientAll;
    }

    @NotNull
    private static Set<String> extractTumorSampleIdsForPatient(@NotNull final String patientIdentifier,
            @NotNull final List<RunContext> runContexts) {
        final Set<String> sampleIdsForPatient = Sets.newHashSet();
        runContexts.forEach(runContext -> {
            final String sampleId = runContext.tumorSample();
            if (sampleId.startsWith(patientIdentifier)) {
                sampleIdsForPatient.add(sampleId);
            }
        });
        return sampleIdsForPatient;
    }

    @NotNull
    private static EcrfModels loadEcrfModels(@NotNull CommandLine cmd) throws IOException, XMLStreamException {
        final String cpctEcrfFilePath = cmd.getOptionValue(CPCT_ECRF_FILE);
        final String cpctFormStatusCsv = cmd.getOptionValue(CPCT_FORM_STATUS_CSV);
        LOGGER.info(String.format("Loading CPCT eCRF from %s.", cpctEcrfFilePath));
        final FormStatusModel cpctFormStatusModel = FormStatusReader.buildModelFromCsv(cpctFormStatusCsv);
        final EcrfModel cpctEcrfModel = EcrfModel.loadFromXMLWithFormStates(cpctEcrfFilePath, cpctFormStatusModel);
        LOGGER.info(String.format("Finished loading CPCT eCRF. Read %s patients.", cpctEcrfModel.patientCount()));

        final String drupEcrfFilePath = cmd.getOptionValue(DRUP_ECRF_FILE);
        LOGGER.info(String.format("Loading DRUP eCRF from %s.", drupEcrfFilePath));
        final EcrfModel drupEcrfModel = EcrfModel.loadFromXMLNoFormStates(drupEcrfFilePath);
        LOGGER.info(String.format("Finished loading DRUP eCRF. Read %s patients.", drupEcrfModel.patientCount()));

        return ImmutableEcrfModels.of(cpctEcrfModel, drupEcrfModel);
    }

    @NotNull
    private static DatabaseAccess createDbWriter(@NotNull CommandLine cmd) throws SQLException {
        final String jdbcUrl = "jdbc:" + cmd.getOptionValue(DB_URL);
        return new DatabaseAccess(cmd.getOptionValue(DB_USER), cmd.getOptionValue(DB_PASS), jdbcUrl);
    }

    private static boolean checkInputs(@NotNull CommandLine cmd) {
        final String runsFolderPathDb = cmd.getOptionValue(RUNS_DIR_DATABASE);
        final String runsFolderPathNonDb = cmd.getOptionValue(RUNS_DIR_NON_DATABASE);

        boolean allParamsPresent = !Utils.anyNull(runsFolderPathDb,
                runsFolderPathNonDb,
                cmd.getOptionValue(DB_USER),
                cmd.getOptionValue(DB_PASS),
                cmd.getOptionValue(DB_URL),
                cmd.getOptionValue(CPCT_ECRF_FILE),
                cmd.getOptionValue(CPCT_FORM_STATUS_CSV),
                cmd.getOptionValue(DRUP_ECRF_FILE),
                cmd.getOptionValue(LIMS_DIRECTORY),
                cmd.getOptionValue(CSV_OUT_DIR));

        boolean validRunDirectories = true;
        if (allParamsPresent) {
            final File runDirectoryDb = new File(runsFolderPathDb);

            if (!runDirectoryDb.exists() || !runDirectoryDb.isDirectory()) {
                validRunDirectories = false;
                LOGGER.warn("HMF database run directory " + runDirectoryDb + " does not exist or is not a directory.");
            }

            final File runDirectoryNonDb = new File(runsFolderPathNonDb);
            if (!runDirectoryNonDb.exists() || !runDirectoryNonDb.isDirectory()) {
                validRunDirectories = false;
                LOGGER.warn("Non-database run directory " + runDirectoryDb + " does not exist or is not a directory.");
            }
        }

        return validRunDirectories && allParamsPresent;
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIR_DATABASE,
                true,
                "Path towards the folder containing patient runs that are considered part of HMF database.");

        options.addOption(RUNS_DIR_NON_DATABASE,
                true,
                "Path towards the folder containing patient runs that are not considered part of HMF database.");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        options.addOption(CPCT_ECRF_FILE, true, "Path towards the cpct ecrf file.");
        options.addOption(CPCT_FORM_STATUS_CSV, true, "Path towards the cpct form status csv file.");
        options.addOption(DRUP_ECRF_FILE, true, "Path towards the drup ecrf file.");
        options.addOption(DO_LOAD_RAW_ECRF, false, "Also write raw ecrf data to database?");

        options.addOption(LIMS_DIRECTORY, true, "Path towards the LIMS directory.");

        options.addOption(CSV_OUT_DIR, true, "Path towards the output directory for csv data dumps.");
        options.addOption(TUMOR_LOCATION_SYMLINK, true, "Name of cancer type csv symlink.");
        options.addOption(PORTAL_DATA_LINK, true, "Name of portal data csv symlink.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
