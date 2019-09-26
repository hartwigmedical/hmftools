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
import org.jetbrains.annotations.Nullable;

public final class LoadClinicalData {

    private static final Logger LOGGER = LogManager.getLogger(LoadClinicalData.class);
    private static final String VERSION = LoadClinicalData.class.getPackage().getImplementationVersion();

    private static final String RUNS_DIRECTORY = "runs_dir";
    private static final String CPCT_ECRF_FILE = "cpct_ecrf";
    private static final String CPCT_FORM_STATUS_CSV = "cpct_form_status_csv";
    private static final String DRUP_ECRF_FILE = "drup_ecrf";
    private static final String DO_LOAD_RAW_ECRF = "do_load_raw_ecrf";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String LIMS_DIRECTORY = "lims_dir";

    private static final String TUMOR_LOCATION_OUTPUT_DIRECTORY = "tumor_location_dir";
    private static final String TUMOR_LOCATION_SYMLINK = "tumor_location_symlink";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, XMLStreamException, SQLException {
        LOGGER.info("Running patient-db v{}", VERSION);
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (checkInputs(cmd)) {
            LOGGER.info("Connecting to database " + cmd.getOptionValue(DB_URL));
            final DatabaseAccess dbWriter = createDbWriter(cmd);

            LOGGER.info("Loading sequence runs from file system.");
            final List<RunContext> runContexts = loadRunContexts(cmd);
            final Map<String, List<String>> sequencedSamplesPerPatient = extractSequencedSamplesFromRunContexts(runContexts);
            final Set<String> sequencedPatientIds = sequencedSamplesPerPatient.keySet();
            final Set<String> sequencedSampleIds = toUniqueSampleIds(sequencedSamplesPerPatient);

            LOGGER.info(String.format(" Loaded sequence runs for %s patient IDs (%s samples).",
                    sequencedPatientIds.size(),
                    sequencedSampleIds.size()));

            LOGGER.info("Loading sample data from LIMS.");
            final Lims lims = LimsFactory.fromLimsDirectory(cmd.getOptionValue(LIMS_DIRECTORY));
            final Map<String, List<SampleData>> limsSampleDataPerPatient = extractAllSamplesFromLims(lims, sequencedSampleIds);
            LOGGER.info(String.format(" Loaded samples for %s patient IDs (%s samples).",
                    limsSampleDataPerPatient.keySet().size(),
                    countValues(limsSampleDataPerPatient)));

            final EcrfModels ecrfModels = loadEcrfModels(cmd);

            if (cmd.hasOption(DO_LOAD_RAW_ECRF)) {
                writeRawEcrf(dbWriter, sequencedPatientIds, ecrfModels);
            }

            writeClinicalData(dbWriter,
                    sequencedPatientIds,
                    limsSampleDataPerPatient,
                    ecrfModels,
                    cmd.getOptionValue(TUMOR_LOCATION_OUTPUT_DIRECTORY),
                    Optional.ofNullable(cmd.getOptionValue(TUMOR_LOCATION_SYMLINK)));
        } else {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db", options);
        }
    }

    @NotNull
    private static List<RunContext> loadRunContexts(@NotNull CommandLine cmd) throws IOException {
        final String runsDirectory = cmd.getOptionValue(RUNS_DIRECTORY);
        final List<RunContext> runContexts = RunsFolderReader.extractRunContexts(new File(runsDirectory));
        LOGGER.info(String.format(" Loaded run contexts from %s (%s sets).", runsDirectory, runContexts.size()));

        return runContexts;
    }

    @NotNull
    private static Map<String, List<String>> extractSequencedSamplesFromRunContexts(@NotNull List<RunContext> runContexts) {
        Map<String, List<String>> sequencedSamplesPerPatient = Maps.newHashMap();
        for (RunContext runContext : runContexts) {
            String patientId = Utils.extractPatientIdentifier(runContext.setName());
            List<String> currentSampleIds = sequencedSamplesPerPatient.get(patientId);
            if (currentSampleIds == null) {
                currentSampleIds = Lists.newArrayList(runContext.tumorSample());
            } else {
                currentSampleIds.add(runContext.tumorSample());
            }
            sequencedSamplesPerPatient.put(patientId, currentSampleIds);
        }

        return sequencedSamplesPerPatient;
    }

    @NotNull
    private static Map<String, List<SampleData>> extractAllSamplesFromLims(@NotNull Lims lims, @NotNull Set<String> sequencedSampleIds) {
        LimsSampleReader sampleReader = new LimsSampleReader(lims, sequencedSampleIds);

        Map<String, List<SampleData>> samplesPerPatient = Maps.newHashMap();
        for (String sampleBarcode : lims.sampleBarcodes()) {
            String sampleId = lims.sampleId(sampleBarcode);
            LimsSampleType sampleType = LimsSampleType.fromSampleId(sampleId);

            if (sampleType != LimsSampleType.OTHER) {
                String patientId = lims.patientId(sampleBarcode);
                SampleData sampleData = sampleReader.read(sampleBarcode, sampleId);

                if (sampleData != null) {
                    List<SampleData> currentSamples = samplesPerPatient.get(patientId);
                    if (currentSamples == null) {
                        currentSamples = Lists.newArrayList(sampleData);
                    } else if (!sampleIdExistsInSampleDataList(currentSamples, sampleId)) {
                        // TODO If a single sample exists in LIMS with multiple barcodes we should pick "the most relevant one".
                        // Currently just picking a random one - sampleId has to be unique in this list.
                        currentSamples.add(sampleData);
                    }
                    samplesPerPatient.put(patientId, currentSamples);
                }
            }
        }

        return samplesPerPatient;
    }

    private static boolean sampleIdExistsInSampleDataList(@NotNull List<SampleData> samples, @NotNull String sampleId) {
        for (SampleData sample : samples) {
            if (sample.sampleId().equals(sampleId)) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private static EcrfModels loadEcrfModels(@NotNull CommandLine cmd) throws IOException, XMLStreamException {
        final String cpctEcrfFilePath = cmd.getOptionValue(CPCT_ECRF_FILE);
        final String cpctFormStatusCsv = cmd.getOptionValue(CPCT_FORM_STATUS_CSV);
        LOGGER.info(String.format("Loading CPCT eCRF from %s.", cpctEcrfFilePath));
        final FormStatusModel cpctFormStatusModel = FormStatusReader.buildModelFromCsv(cpctFormStatusCsv);
        final EcrfModel cpctEcrfModel = EcrfModel.loadFromXMLWithFormStates(cpctEcrfFilePath, cpctFormStatusModel);
        LOGGER.info(String.format(" Finished loading CPCT eCRF. Read %s patients.", cpctEcrfModel.patientCount()));

        final String drupEcrfFilePath = cmd.getOptionValue(DRUP_ECRF_FILE);
        LOGGER.info(String.format("Loading DRUP eCRF from %s.", drupEcrfFilePath));
        final EcrfModel drupEcrfModel = EcrfModel.loadFromXMLNoFormStates(drupEcrfFilePath);
        LOGGER.info(String.format(" Finished loading DRUP eCRF. Read %s patients.", drupEcrfModel.patientCount()));

        return ImmutableEcrfModels.of(cpctEcrfModel, drupEcrfModel);
    }

    private static void writeRawEcrf(@NotNull DatabaseAccess dbWriter, @NotNull Set<String> sequencedPatients,
            @NotNull EcrfModels ecrfModels) {
        final EcrfModel cpctEcrfModel = ecrfModels.cpctModel();
        LOGGER.info(String.format("Writing raw cpct ecrf data for %s patients", cpctEcrfModel.patientCount()));
        dbWriter.clearCpctEcrf();
        dbWriter.writeCpctEcrf(cpctEcrfModel, sequencedPatients);
        LOGGER.info(String.format(" Finished writing raw cpct ecrf data for %s patients.", cpctEcrfModel.patientCount()));

        final EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info(String.format("Writing raw drup ecrf data for %s patients", drupEcrfModel.patientCount()));
        dbWriter.clearDrupEcrf();
        dbWriter.writeDrupEcrf(drupEcrfModel, sequencedPatients);
        LOGGER.info(String.format(" Finished writing raw drup ecrf data for %s patients.", drupEcrfModel.patientCount()));
    }

    private static void writeClinicalData(@NotNull DatabaseAccess dbAccess, @NotNull Set<String> sequencedPatientIds,
            @NotNull Map<String, List<SampleData>> limsSampleDataPerPatient, @NotNull EcrfModels ecrfModels,
            @NotNull String tumorLocationOutputDir, @NotNull Optional<String> tumorLocationSymlink) throws IOException {
        TumorLocationCurator tumorLocationCurator = TumorLocationCurator.fromProductionResource();
        BiopsySiteCurator biopsySiteCurator = BiopsySiteCurator.fromProductionResource();
        TreatmentCurator treatmentCurator = TreatmentCurator.fromProductionResource();

        Map<String, Patient> patients =
                loadAndInterpretPatients(limsSampleDataPerPatient, ecrfModels, tumorLocationCurator, treatmentCurator, biopsySiteCurator);

        DumpTumorLocationData.writeCuratedTumorLocationsToCSV(tumorLocationOutputDir, tumorLocationSymlink, patients.values());

        LOGGER.info("Clearing interpreted clinical tables in database.");
        dbAccess.clearClinicalTables();

        int missingPatients = 0;
        int missingSamples = 0;
        LOGGER.info(String.format("Writing clinical data for %s sequenced patients.", sequencedPatientIds.size()));
        for (final String patientId : sequencedPatientIds) {
            Patient patient = patients.get(patientId);
            if (patient == null) {
                missingPatients++;
                List<SampleData> samples = limsSampleDataPerPatient.get(patientId);
                if (samples == null) {
                    LOGGER.warn("Could not find any samples for " + patientId + "! Skipping writing to db.");
                } else {
                    List<SampleData> sequencedSamples = sequencedOnly(samples);
                    missingSamples += sequencedSamples.size();
                    dbAccess.writeSampleClinicalData(patientId, sequencedSamples);
                }
            } else if (patient.sequencedBiopsies().isEmpty()) {
                LOGGER.warn("No sequenced biopsies found for sequenced patient: " + patientId + "! Skipping writing to db.");
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
    private static Map<String, Patient> loadAndInterpretPatients(@NotNull Map<String, List<SampleData>> limsSampleDataPerPatient,
            @NotNull EcrfModels ecrfModels, @NotNull TumorLocationCurator tumorLocationCurator, @NotNull TreatmentCurator treatmentCurator,
            @NotNull BiopsySiteCurator biopsySiteCurator) {
        final EcrfModel cpctEcrfModel = ecrfModels.cpctModel();
        LOGGER.info(String.format("Interpreting and curating data for %s CPCT patients.", cpctEcrfModel.patientCount()));
        EcrfPatientReader cpctPatientReader = new CpctPatientReader(tumorLocationCurator,
                CpctUtil.extractHospitalMap(cpctEcrfModel),
                biopsySiteCurator,
                treatmentCurator);

        Map<String, Patient> cpctPatients = readEcrfPatients(cpctPatientReader, cpctEcrfModel.patients(), limsSampleDataPerPatient);
        LOGGER.info(String.format(" Finished curation of %s CPCT patients.", cpctPatients.size()));

        final EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info(String.format("Interpreting and curating data for %s DRUP patients.", drupEcrfModel.patientCount()));
        EcrfPatientReader drupPatientReader = new DrupPatientReader(tumorLocationCurator, biopsySiteCurator);

        Map<String, Patient> drupPatients = readEcrfPatients(drupPatientReader, drupEcrfModel.patients(), limsSampleDataPerPatient);
        LOGGER.info(String.format(" Finished curation of %s DRUP patients.", drupPatients.size()));

        LOGGER.info("Interpreting and curating data based off LIMS (WIDE and CORE).");
        Map<String, Patient> patientsFromLims = readLimsPatients(limsSampleDataPerPatient, tumorLocationCurator);
        LOGGER.info(String.format(" Finished curation of %s patients based off LIMS", patientsFromLims.keySet().size()));

        Map<String, Patient> mergedPatients = Maps.newHashMap();
        mergedPatients.putAll(cpctPatients);
        mergedPatients.putAll(drupPatients);
        mergedPatients.putAll(patientsFromLims);
        mergedPatients.putAll(readColoPatients());
        return mergedPatients;
    }

    @NotNull
    private static Map<String, Patient> readEcrfPatients(@NotNull EcrfPatientReader reader, @NotNull Iterable<EcrfPatient> patients,
            @NotNull Map<String, List<SampleData>> limsSampleDataPerPatient) {
        final Map<String, Patient> patientMap = Maps.newHashMap();
        for (final EcrfPatient ecrfPatient : patients) {
            List<SampleData> sequencedSamples = sequencedOnly(limsSampleDataPerPatient.get(ecrfPatient.patientId()));
            Patient patient = reader.read(ecrfPatient, sequencedSamples);
            patientMap.put(patient.patientIdentifier(), patient);
        }
        return patientMap;
    }

    @NotNull
    private static Map<String, Patient> readLimsPatients(@NotNull Map<String, List<SampleData>> limsSampleDataPerPatient,
            @NotNull TumorLocationCurator tumorLocationCurator) {
        final Map<String, Patient> patientMap = Maps.newHashMap();
        final LimsPatientReader limsPatientReader = new LimsPatientReader(tumorLocationCurator);

        for (Map.Entry<String, List<SampleData>> entry : limsSampleDataPerPatient.entrySet()) {
            List<SampleData> samples = entry.getValue();

            assert samples != null;
            List<SampleData> tumorSamples = extractTumorSamples(samples);
            if (!tumorSamples.isEmpty()) {
                LimsSampleType sampleType = LimsSampleType.fromSampleId(tumorSamples.get(0).sampleId());

                if (sampleType == LimsSampleType.CORE || sampleType == LimsSampleType.WIDE) {
                    String patientId = entry.getKey();
                    Patient limsPatient =
                            limsPatientReader.read(patientId, tumorSamples.get(0).limsPrimaryTumor(), sequencedOnly(tumorSamples));
                    patientMap.put(patientId, limsPatient);
                }
            }
        }

        return patientMap;
    }

    @NotNull
    private static List<SampleData> extractTumorSamples(@NotNull Iterable<SampleData> samples) {
        List<SampleData> tumorSamples = Lists.newArrayList();

        for (SampleData sample : samples) {
            LimsSampleType sampleType = LimsSampleType.fromSampleId(sample.sampleId());
            if (sampleType != LimsSampleType.OTHER) {
                if (sample.sampleId().substring(12).contains("T")) {
                    tumorSamples.add(sample);
                }
            }
        }
        return tumorSamples;
    }

    @NotNull
    private static Map<String, Patient> readColoPatients() {
        final Map<String, Patient> patientMap = Maps.newHashMap();
        final ColoPatientReader coloPatientReader = new ColoPatientReader();
        LOGGER.info("Creating patient representation for COLO829.");
        Patient colo829Patient = coloPatientReader.read("COLO829T");

        patientMap.put(colo829Patient.patientIdentifier(), colo829Patient);
        return patientMap;
    }

    private static <V, K> int countValues(@NotNull Map<V, List<K>> map) {
        int count = 0;
        for (Map.Entry<V, List<K>> entry : map.entrySet()) {
            count += entry.getValue().size();
        }
        return count;
    }

    @NotNull
    private static Set<String> toUniqueSampleIds(@NotNull Map<String, List<String>> samplesPerPatient) {
        Set<String> uniqueSampleIds = Sets.newHashSet();
        for (Map.Entry<String, List<String>> entry : samplesPerPatient.entrySet()) {
            uniqueSampleIds.addAll(entry.getValue());
        }
        return uniqueSampleIds;
    }

    @NotNull
    private static List<SampleData> sequencedOnly(@Nullable Iterable<SampleData> samples) {
        if (samples == null) {
            return Lists.newArrayList();
        }

        List<SampleData> sequencedSamples = Lists.newArrayList();
        for (SampleData sample : samples) {
            if (sample.sequenced()) {
                sequencedSamples.add(sample);
            }
        }
        return sequencedSamples;
    }

    @NotNull
    private static DatabaseAccess createDbWriter(@NotNull CommandLine cmd) throws SQLException {
        final String jdbcUrl = "jdbc:" + cmd.getOptionValue(DB_URL);
        return new DatabaseAccess(cmd.getOptionValue(DB_USER), cmd.getOptionValue(DB_PASS), jdbcUrl);
    }

    private static boolean checkInputs(@NotNull CommandLine cmd) {
        final String runsDirectory = cmd.getOptionValue(RUNS_DIRECTORY);

        boolean allParamsPresent = !Utils.anyNull(runsDirectory,
                cmd.getOptionValue(DB_USER),
                cmd.getOptionValue(DB_PASS),
                cmd.getOptionValue(DB_URL),
                cmd.getOptionValue(CPCT_ECRF_FILE),
                cmd.getOptionValue(CPCT_FORM_STATUS_CSV),
                cmd.getOptionValue(DRUP_ECRF_FILE),
                cmd.getOptionValue(LIMS_DIRECTORY),
                cmd.getOptionValue(TUMOR_LOCATION_OUTPUT_DIRECTORY));

        boolean validRunDirectories = true;
        if (allParamsPresent) {
            final File runDirectoryDb = new File(runsDirectory);

            if (!runDirectoryDb.exists() || !runDirectoryDb.isDirectory()) {
                validRunDirectories = false;
                LOGGER.warn("HMF database run directory " + runDirectoryDb + " does not exist or is not a directory.");
            }
        }

        return validRunDirectories && allParamsPresent;
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIRECTORY,
                true,
                "Path towards the folder containing patient runs that are considered part of HMF database.");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        options.addOption(CPCT_ECRF_FILE, true, "Path towards the cpct ecrf file.");
        options.addOption(CPCT_FORM_STATUS_CSV, true, "Path towards the cpct form status csv file.");
        options.addOption(DRUP_ECRF_FILE, true, "Path towards the drup ecrf file.");
        options.addOption(DO_LOAD_RAW_ECRF, false, "Also write raw ecrf data to database?");

        options.addOption(LIMS_DIRECTORY, true, "Path towards the LIMS directory.");

        options.addOption(TUMOR_LOCATION_OUTPUT_DIRECTORY, true, "Path towards the output directory for tumor location data dumps.");
        options.addOption(TUMOR_LOCATION_SYMLINK, true, "Name of tumor location csv symlink.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
