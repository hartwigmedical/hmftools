package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_URL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.clinical.PatientTumorCurationStatus;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusModel;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusReader;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.reportingdb.ReportingDatabase;
import com.hartwig.hmftools.common.reportingdb.ReportingEntry;
import com.hartwig.hmftools.patientdb.context.RunContext;
import com.hartwig.hmftools.patientdb.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.readers.ColoPatientReader;
import com.hartwig.hmftools.patientdb.readers.CorePatientReader;
import com.hartwig.hmftools.patientdb.readers.EcrfPatientReader;
import com.hartwig.hmftools.patientdb.readers.LimsSampleReader;
import com.hartwig.hmftools.patientdb.readers.RunsFolderReader;
import com.hartwig.hmftools.patientdb.readers.WidePatientReader;
import com.hartwig.hmftools.patientdb.readers.cpct.CpctPatientReader;
import com.hartwig.hmftools.patientdb.readers.cpct.CpctUtil;
import com.hartwig.hmftools.patientdb.readers.drup.DrupPatientReader;
import com.hartwig.hmftools.patientdb.readers.wide.ImmutableWideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WideAvlTreatmentData;
import com.hartwig.hmftools.patientdb.readers.wide.WideBiopsyData;
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfFileReader;
import com.hartwig.hmftools.patientdb.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.readers.wide.WideFiveDays;
import com.hartwig.hmftools.patientdb.readers.wide.WidePreAvlTreatmentData;
import com.hartwig.hmftools.patientdb.readers.wide.WideResponseData;
import com.hartwig.hmftools.patientdb.validators.CurationValidator;
import com.hartwig.hmftools.patientdb.validators.PatientValidator;

import org.apache.commons.cli.CommandLine;
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
    private static final String PIPELINE_VERSION = "pipeline_version_file";

    private static final String RUNS_DIRECTORY = "runs_dir";

    private static final String CPCT_ECRF_FILE = "cpct_ecrf";
    private static final String CPCT_FORM_STATUS_CSV = "cpct_form_status_csv";
    private static final String DRUP_ECRF_FILE = "drup_ecrf";

    private static final String DO_LOAD_CLINICAL_DATA = "do_load_clinical_data";
    private static final String DO_LOAD_RAW_ECRF = "do_load_raw_ecrf";

    private static final String LIMS_DIRECTORY = "lims_dir";
    private static final String REPORTING_DB_TSV = "reporting_db_tsv";

    private static final String CURATED_PRIMARY_TUMOR_TSV = "curated_primary_tumor_tsv";
    private static final String PATIENT_TUMOR_CURATION_STATUS_TSV = "patient_tumor_curation_status_tsv";

    private static final String DO_PROCESS_WIDE_CLINICAL_DATA = "do_process_wide_clinical_data";
    private static final String WIDE_PRE_AVL_TREATMENT_CSV = "wide_pre_avl_treatment_csv";
    private static final String WIDE_BIOPSY_CSV = "wide_biopsy_csv";
    private static final String WIDE_AVL_TREATMENT_CSV = "wide_avl_treatment_csv";
    private static final String WIDE_RESPONSE_CSV = "wide_response_csv";
    private static final String WIDE_FIVE_DAYS_CSV = "wide_five_days_csv";

    private static final String DOID_JSON = "doid_json";
    private static final String TUMOR_LOCATION_MAPPING_TSV = "tumor_location_mapping_tsv";
    private static final String TREATMENT_MAPPING_CSV = "treatment_mapping_csv";
    private static final String BIOPSY_MAPPING_CSV = "biopsy_mapping_csv";

    public static void main(@NotNull String[] args) throws ParseException, IOException, XMLStreamException, SQLException {
        LOGGER.info("Running Clinical Patient DB v{}", VERSION);
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        if (!checkInputs(cmd)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db", options);
            System.exit(1);
        }

        List<DoidNode> doidNodes = DiseaseOntology.readDoidOwlEntryFromDoidJson(cmd.getOptionValue(DOID_JSON)).nodes();
        PrimaryTumorCurator primaryTumorCurator = new PrimaryTumorCurator(cmd.getOptionValue(TUMOR_LOCATION_MAPPING_TSV), doidNodes);
        BiopsySiteCurator biopsySiteCurator = new BiopsySiteCurator(cmd.getOptionValue(BIOPSY_MAPPING_CSV));
        TreatmentCurator treatmentCurator = new TreatmentCurator(cmd.getOptionValue(TREATMENT_MAPPING_CSV));

        LOGGER.info("Loading sequence runs from {}", cmd.getOptionValue(RUNS_DIRECTORY));
        List<RunContext> runContexts = loadRunContexts(cmd.getOptionValue(RUNS_DIRECTORY), cmd.getOptionValue(PIPELINE_VERSION));
        Map<String, List<String>> sequencedSamplesPerPatient = extractSequencedSamplesFromRunContexts(runContexts);
        Map<String, String> sampleToSetNameMap = extractSampleToSetNameMap(runContexts);
        Set<String> sequencedPatientIds = sequencedSamplesPerPatient.keySet();

        LOGGER.info(" Loaded sequence runs for {} patient IDs ({} samples)",
                sequencedPatientIds.size(),
                toUniqueSampleIds(sequencedSamplesPerPatient).size());

        LOGGER.info("Loading sample data from LIMS in {}", cmd.getOptionValue(LIMS_DIRECTORY));
        Lims lims = LimsFactory.fromLimsDirectory(cmd.getOptionValue(LIMS_DIRECTORY));
        Map<String, List<SampleData>> sampleDataPerPatient =
                extractAllSamplesFromLims(lims, sampleToSetNameMap, sequencedSamplesPerPatient);
        LOGGER.info(" Loaded samples for {} patient IDs ({} samples)",
                sampleDataPerPatient.keySet().size(),
                countValues(sampleDataPerPatient));

        EcrfModels ecrfModels = loadEcrfModels(cmd);

        List<Patient> patients = interpret(sampleDataPerPatient, ecrfModels, primaryTumorCurator, biopsySiteCurator, treatmentCurator);

        LOGGER.info("Check for missing curation tumor location when info is known");
        Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap =
                generatePatientTumorCurationStatusMap(sampleDataPerPatient, patients, cmd.getOptionValue(REPORTING_DB_TSV));

        LOGGER.info("Writing patient tumor curation status");
        DumpPrimaryTumorData.writePatientTumorCurationStatesToTSV(cmd.getOptionValue(PATIENT_TUMOR_CURATION_STATUS_TSV),
                patientTumorCurationStatusMap);

        LOGGER.info("Writing curated primary tumors");
        DumpPrimaryTumorData.writeCuratedPrimaryTumorsToTSV(cmd.getOptionValue(CURATED_PRIMARY_TUMOR_TSV), patients);

        if (cmd.hasOption(DO_LOAD_CLINICAL_DATA)) {
            LOGGER.info("Connecting to database {}", cmd.getOptionValue(DB_URL));
            DatabaseAccess dbWriter = databaseAccess(cmd);

            if (cmd.hasOption(DO_LOAD_RAW_ECRF)) {
                writeRawEcrf(dbWriter, sequencedPatientIds, ecrfModels);
            }

            writeClinicalData(dbWriter, lims, sequencedPatientIds, sampleDataPerPatient, patients);

            dbWriter.writeValidationFindings(CurationValidator.validatePrimaryTumorCurator(primaryTumorCurator));
            dbWriter.writeValidationFindings(CurationValidator.validateTreatmentCurator(treatmentCurator));
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Map<String, PatientTumorCurationStatus> generatePatientTumorCurationStatusMap(
            @NotNull Map<String, List<SampleData>> sampleDataPerPatient, @NotNull List<Patient> patients, @NotNull String reportingDbTsv)
            throws IOException {
        Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap = Maps.newHashMap();

        List<String> reportedBarcodes = Lists.newArrayList();
        for (ReportingEntry entry : ReportingDatabase.read(reportingDbTsv)) {
            reportedBarcodes.add(entry.tumorBarcode());
        }

        List<String> patientsWithSamplesToBeReported = Lists.newArrayList();
        for (Map.Entry<String, List<SampleData>> sampleDataEntry : sampleDataPerPatient.entrySet()) {
            String patientId = sampleDataEntry.getKey();
            for (SampleData sampleData : sampleDataEntry.getValue()) {
                if (!sampleData.requiresCuratedPrimaryTumor()) {
                    patientTumorCurationStatusMap.put(patientId, PatientTumorCurationStatus.NEEDS_NO_CURATED_PRIMARY_TUMOR);
                } else if (sampleData.isSomaticTumorSample()) {
                    if (reportedBarcodes.contains(sampleData.sampleBarcode())) {
                        if (!patientId.startsWith("CORE") && !patientId.startsWith("WIDE")) {
                            patientTumorCurationStatusMap.put(patientId, PatientTumorCurationStatus.ALREADY_REPORTED);
                        }
                    } else {
                        patientsWithSamplesToBeReported.add(patientId);
                    }
                }
            }
        }

        for (String patientId : patientsWithSamplesToBeReported) {
            Patient patient = findByPatientId(patients, patientId);
            if (patient != null) {
                String tumorLocationSearchTerm = patient.baselineData().curatedPrimaryTumor().searchTerm();
                if (tumorLocationSearchTerm != null && !tumorLocationSearchTerm.isEmpty()) {
                    if (patient.baselineData().curatedPrimaryTumor().location() == null) {
                        LOGGER.warn("Could not curate patient {} for primary tumor '{}'",
                                patient.patientIdentifier(),
                                tumorLocationSearchTerm);
                    }
                } else {
                    if (patient.patientIdentifier().startsWith("CORE") || patient.patientIdentifier().startsWith("WIDE")) {
                        LOGGER.warn("Could not find input tumor location for patient {}", patient.patientIdentifier());
                    } else {
                        patientTumorCurationStatusMap.put(patient.patientIdentifier(), PatientTumorCurationStatus.MISSING_TUMOR_CURATION);
                    }
                }
            } else {
                patientTumorCurationStatusMap.put(patientId, PatientTumorCurationStatus.NOT_RESOLVED);
            }
        }

        return patientTumorCurationStatusMap;
    }

    @NotNull
    private static List<RunContext> loadRunContexts(@NotNull String runsDirectory, @NotNull String pipelineVersionFile) throws IOException {
        List<RunContext> runContexts = RunsFolderReader.extractRunContexts(new File(runsDirectory), pipelineVersionFile);
        LOGGER.info(" Loaded run contexts from {} ({} sets)", runsDirectory, runContexts.size());

        return runContexts;
    }

    @NotNull
    private static Map<String, String> extractSampleToSetNameMap(@NotNull List<RunContext> runContexts) {
        Map<String, String> sampleToSetNameMap = Maps.newHashMap();
        for (RunContext runContext : runContexts) {
            if (sampleToSetNameMap.containsKey(runContext.tumorSample())) {
                LOGGER.warn("Duplicate sample ID found in run contexts: {}", runContext.tumorSample());
            }
            sampleToSetNameMap.put(runContext.tumorSample(), runContext.setName());
        }
        return sampleToSetNameMap;
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
    private static Map<String, List<SampleData>> extractAllSamplesFromLims(@NotNull Lims lims,
            @NotNull Map<String, String> sampleToSetNameMap, @NotNull Map<String, List<String>> sequencedSamplesPerPatient) {
        LimsSampleReader sampleReader = new LimsSampleReader(lims, sampleToSetNameMap, toUniqueSampleIds(sequencedSamplesPerPatient));

        Map<String, List<SampleData>> samplesPerPatient = Maps.newHashMap();
        for (String sampleBarcode : lims.sampleBarcodes()) {
            String sampleId = lims.sampleId(sampleBarcode);
            SampleData sampleData = sampleReader.read(sampleBarcode, sampleId);

            if (sampleData != null) {
                String patientId = lims.patientId(sampleBarcode);
                List<SampleData> currentSamples = samplesPerPatient.get(patientId);
                if (currentSamples == null) {
                    currentSamples = Lists.newArrayList(sampleData);
                } else if (!sampleIdExistsInSampleDataList(currentSamples, sampleId)) {
                    // Ideally if a single sample exists in LIMS with multiple barcodes we pick "the most relevant one".
                    // Currently just picking a random one - sampleId has to be unique in this list.
                    currentSamples.add(sampleData);
                }
                samplesPerPatient.put(patientId, currentSamples);
            }
        }

        // Some samples may be missing from LIMS simply because they are old and we did not collect information back in those days.
        for (Map.Entry<String, List<String>> sequencedPatientEntry : sequencedSamplesPerPatient.entrySet()) {
            List<SampleData> samples = samplesPerPatient.get(sequencedPatientEntry.getKey());
            if (samples == null) {
                samples = Lists.newArrayList();
            }
            for (String sampleId : sequencedPatientEntry.getValue()) {
                if (!sampleIdExistsInSampleDataList(samples, sampleId)) {
                    LOGGER.info(" Creating sample data for {}. This sample is not found in LIMS even though it has been sequenced!",
                            sampleId);
                    SampleData sampleData = sampleReader.readSequencedSampleWithoutBarcode(sampleId);
                    if (sampleData != null) {
                        samples.add(sampleData);
                    }
                }
            }
            samplesPerPatient.put(sequencedPatientEntry.getKey(), samples);
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
        EcrfModel cpctEcrfModel = buildCpctEcrfModel(cmd);
        EcrfModel drupEcrfModel = buildDrupEcrfModel(cmd);
        WideEcrfModel wideEcrfModel = buildWideEcrfModel(cmd);

        return ImmutableEcrfModels.builder().cpctModel(cpctEcrfModel).drupModel(drupEcrfModel).wideModel(wideEcrfModel).build();
    }

    @NotNull
    private static EcrfModel buildCpctEcrfModel(@NotNull CommandLine cmd) throws IOException, XMLStreamException {
        String cpctEcrfFilePath = cmd.getOptionValue(CPCT_ECRF_FILE);
        String cpctFormStatusCsv = cmd.getOptionValue(CPCT_FORM_STATUS_CSV);
        LOGGER.info("Loading CPCT eCRF from {}", cpctEcrfFilePath);
        FormStatusModel cpctFormStatusModel = FormStatusReader.buildModelFromCsv(cpctFormStatusCsv);
        EcrfModel cpctEcrfModel = EcrfModel.loadFromXMLWithFormStates(cpctEcrfFilePath, cpctFormStatusModel);
        LOGGER.info(" Finished loading CPCT eCRF. Read {} patients", cpctEcrfModel.patientCount());

        return cpctEcrfModel;
    }

    @NotNull
    private static EcrfModel buildDrupEcrfModel(@NotNull CommandLine cmd) throws FileNotFoundException, XMLStreamException {
        String drupEcrfFilePath = cmd.getOptionValue(DRUP_ECRF_FILE);
        LOGGER.info("Loading DRUP eCRF from {}", drupEcrfFilePath);
        EcrfModel drupEcrfModel = EcrfModel.loadFromXMLNoFormStates(drupEcrfFilePath);
        LOGGER.info(" Finished loading DRUP eCRF. Read {} patients", drupEcrfModel.patientCount());

        return drupEcrfModel;
    }

    @NotNull
    private static WideEcrfModel buildWideEcrfModel(@NotNull CommandLine cmd) throws IOException {
        WideEcrfModel wideEcrfModel;

        if (cmd.hasOption(DO_PROCESS_WIDE_CLINICAL_DATA)) {
            LOGGER.info("Loading WIDE eCRF");

            String preAvlTreatmentCsv = cmd.getOptionValue(WIDE_PRE_AVL_TREATMENT_CSV);
            List<WidePreAvlTreatmentData> preAvlTreatments = WideEcrfFileReader.readPreAvlTreatments(preAvlTreatmentCsv);
            LOGGER.info(" Loaded {} WIDE pre-AVL-treatments from {}", preAvlTreatments.size(), preAvlTreatmentCsv);

            String biopsyCsv = cmd.getOptionValue(WIDE_BIOPSY_CSV);
            List<WideBiopsyData> biopsies = WideEcrfFileReader.readBiopsies(biopsyCsv);
            LOGGER.info(" Loaded {} WIDE biopsies from {}", biopsies.size(), biopsyCsv);

            String avlTreatmentCsv = cmd.getOptionValue(WIDE_AVL_TREATMENT_CSV);
            List<WideAvlTreatmentData> avlTreatments = WideEcrfFileReader.readAvlTreatments(avlTreatmentCsv);
            LOGGER.info(" Loaded {} WIDE AVL treatments from {}", avlTreatments.size(), avlTreatmentCsv);

            String wideResponseCsv = cmd.getOptionValue(WIDE_RESPONSE_CSV);
            List<WideResponseData> responses = WideEcrfFileReader.readResponses(wideResponseCsv);
            LOGGER.info(" Loaded {} WIDE responses from {}", responses.size(), wideResponseCsv);

            String fiveDaysCsv = cmd.getOptionValue(WIDE_FIVE_DAYS_CSV);
            List<WideFiveDays> fiveDays = WideEcrfFileReader.readFiveDays(fiveDaysCsv);
            LOGGER.info(" Loaded {} WIDE five days entries from {}", fiveDays.size(), fiveDaysCsv);

            wideEcrfModel = ImmutableWideEcrfModel.builder()
                    .preAvlTreatments(preAvlTreatments)
                    .biopsies(biopsies)
                    .avlTreatments(avlTreatments)
                    .responses(responses)
                    .fiveDays(fiveDays)
                    .build();
        } else {
            LOGGER.info("Skipping the loading of WIDE eCRF");
            wideEcrfModel = ImmutableWideEcrfModel.builder()
                    .preAvlTreatments(Lists.newArrayList())
                    .biopsies(Lists.newArrayList())
                    .avlTreatments(Lists.newArrayList())
                    .responses(Lists.newArrayList())
                    .fiveDays(Lists.newArrayList())
                    .build();
        }

        return wideEcrfModel;
    }

    private static void writeRawEcrf(@NotNull DatabaseAccess dbWriter, @NotNull Set<String> sequencedPatients,
            @NotNull EcrfModels ecrfModels) {
        EcrfModel cpctEcrfModel = ecrfModels.cpctModel();
        LOGGER.info("Writing raw cpct ecrf data for {} patients", cpctEcrfModel.patientCount());
        dbWriter.clearCpctEcrf();
        dbWriter.writeCpctEcrf(cpctEcrfModel, sequencedPatients);
        LOGGER.info(" Finished writing raw cpct ecrf data for {} patients", cpctEcrfModel.patientCount());

        EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info("Writing raw drup ecrf data for {} patients", drupEcrfModel.patientCount());
        dbWriter.clearDrupEcrf();
        dbWriter.writeDrupEcrf(drupEcrfModel, sequencedPatients);
        LOGGER.info(" Finished writing raw drup ecrf data for {} patients", drupEcrfModel.patientCount());
    }

    @NotNull
    private static List<Patient> interpret(@NotNull Map<String, List<SampleData>> sampleDataPerPatient, @NotNull EcrfModels ecrfModels,
            @NotNull PrimaryTumorCurator primaryTumorCurator, @NotNull BiopsySiteCurator biopsySiteCurator,
            @NotNull TreatmentCurator treatmentCurator) {
        EcrfModel cpctEcrfModel = ecrfModels.cpctModel();
        LOGGER.info("Interpreting and curating data for {} CPCT patients", cpctEcrfModel.patientCount());
        EcrfPatientReader cpctPatientReader =
                new CpctPatientReader(primaryTumorCurator, CpctUtil.extractHospitalMap(cpctEcrfModel), biopsySiteCurator, treatmentCurator);

        List<Patient> cpctPatients = readEcrfPatients(cpctPatientReader, cpctEcrfModel.patients(), sampleDataPerPatient);
        LOGGER.info(" Finished curation of {} CPCT patients", cpctPatients.size());

        EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info("Interpreting and curating data for {} DRUP patients", drupEcrfModel.patientCount());
        EcrfPatientReader drupPatientReader = new DrupPatientReader(primaryTumorCurator, biopsySiteCurator);

        List<Patient> drupPatients = readEcrfPatients(drupPatientReader, drupEcrfModel.patients(), sampleDataPerPatient);
        LOGGER.info(" Finished curation of {} DRUP patients", drupPatients.size());

        LOGGER.info("Interpreting and curating data for WIDE patients");
        List<Patient> widePatients = readWidePatients(ecrfModels.wideModel(), sampleDataPerPatient, primaryTumorCurator, treatmentCurator);
        LOGGER.info(" Finished curation of {} WIDE patients", widePatients.size());

        LOGGER.info("Interpreting and curating data for CORE patients");
        List<Patient> corePatients = readCorePatients(sampleDataPerPatient, primaryTumorCurator);
        LOGGER.info(" Finished curation of {} CORE patients", corePatients.size());

        List<Patient> mergedPatients = Lists.newArrayList();
        mergedPatients.addAll(cpctPatients);
        mergedPatients.addAll(drupPatients);
        mergedPatients.addAll(widePatients);
        mergedPatients.addAll(corePatients);
        mergedPatients.addAll(readColoPatients());
        return mergedPatients;
    }

    @NotNull
    private static List<Patient> readEcrfPatients(@NotNull EcrfPatientReader reader, @NotNull Iterable<EcrfPatient> ecrfPatients,
            @NotNull Map<String, List<SampleData>> sampleDataPerPatient) {
        List<Patient> patients = Lists.newArrayList();
        for (EcrfPatient ecrfPatient : ecrfPatients) {
            List<SampleData> sequencedSamples = sequencedOnly(sampleDataPerPatient.get(ecrfPatient.patientId()));
            patients.add(reader.read(ecrfPatient, sequencedSamples));
        }
        return patients;
    }

    @NotNull
    private static List<Patient> readWidePatients(@NotNull WideEcrfModel wideEcrfModel,
            @NotNull Map<String, List<SampleData>> sampleDataPerPatient, @NotNull PrimaryTumorCurator primaryTumorCurator,
            @NotNull TreatmentCurator treatmentCurator) {
        List<Patient> patients = Lists.newArrayList();

        WidePatientReader widePatientReader = new WidePatientReader(wideEcrfModel, primaryTumorCurator, treatmentCurator);
        for (Map.Entry<String, List<SampleData>> entry : sampleDataPerPatient.entrySet()) {
            List<SampleData> tumorSamples = tumorSamplesOnly(entry.getValue());
            if (!tumorSamples.isEmpty() && tumorSamples.get(0).cohortId().equals("WIDE")) {
                String patientId = entry.getKey();
                // We assume every sample for a single patient has the same primary tumor.
                String primaryTumor = tumorSamples.get(0).limsPrimaryTumor();
                patients.add(widePatientReader.read(patientId, primaryTumor, sequencedOnly(tumorSamples)));
            }
        }
        return patients;
    }

    @NotNull
    private static List<Patient> readCorePatients(@NotNull Map<String, List<SampleData>> sampleDataPerPatient,
            @NotNull PrimaryTumorCurator primaryTumorCurator) {
        List<Patient> patients = Lists.newArrayList();
        CorePatientReader corePatientReader = new CorePatientReader(primaryTumorCurator);

        for (Map.Entry<String, List<SampleData>> entry : sampleDataPerPatient.entrySet()) {
            List<SampleData> tumorSamples = tumorSamplesOnly(entry.getValue());
            if (!tumorSamples.isEmpty() && tumorSamples.get(0).cohortId().contains("CORE")) {
                String patientId = entry.getKey();
                // We assume every sample for a single patient has the same primary tumor.
                String primaryTumor = tumorSamples.get(0).limsPrimaryTumor();
                patients.add(corePatientReader.read(patientId, primaryTumor, sequencedOnly(tumorSamples)));
            }
        }

        return patients;
    }

    @NotNull
    private static List<Patient> readColoPatients() {
        List<Patient> patients = Lists.newArrayList();

        ColoPatientReader coloPatientReader = new ColoPatientReader();
        LOGGER.info("Creating patient representation for COLO829");
        patients.add(coloPatientReader.read("COLO829T"));

        return patients;
    }

    private static void writeClinicalData(@NotNull DatabaseAccess dbAccess, @NotNull Lims lims, @NotNull Set<String> sequencedPatientIds,
            @NotNull Map<String, List<SampleData>> sampleDataPerPatient, @NotNull List<Patient> patients) {
        LOGGER.info("Clearing interpreted clinical tables in database");
        dbAccess.clearClinicalTables();

        int missingPatients = 0;
        int missingSamples = 0;
        LOGGER.info("Writing clinical data for {} sequenced patients", sequencedPatientIds.size());
        for (String patientId : sequencedPatientIds) {
            Patient patient = findByPatientId(patients, patientId);
            if (patient == null) {
                LOGGER.warn("No clinical data found for patient {}", patientId);
                missingPatients++;
                List<SampleData> sequencedSamples = sequencedOnly(sampleDataPerPatient.get(patientId));
                missingSamples += sequencedSamples.size();
                dbAccess.writeSampleClinicalData(patientId, lims.isBlacklisted(patientId), sequencedSamples);
            } else if (patient.sequencedBiopsies().isEmpty()) {
                LOGGER.warn("No sequenced biopsies found for sequenced patient: {}! Skipping writing to db", patientId);
            } else {
                dbAccess.writeFullClinicalData(patient, lims.isBlacklisted(patientId));
                List<ValidationFinding> findings = PatientValidator.validatePatient(patient);

                dbAccess.writeValidationFindings(findings);
                dbAccess.writeValidationFindings(patient.matchFindings());
            }
        }

        if (missingPatients > 0) {
            LOGGER.warn("Could not load {} patients ({} samples)!", missingPatients, missingSamples);
        }
    }

    @Nullable
    private static Patient findByPatientId(@NotNull List<Patient> patients, @NotNull String patientId) {
        for (Patient patient : patients) {
            if (patient.patientIdentifier().equals(patientId)) {
                return patient;
            }
        }

        return null;
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
    private static List<SampleData> tumorSamplesOnly(@NotNull Iterable<SampleData> samples) {
        List<SampleData> tumorSamples = Lists.newArrayList();

        for (SampleData sample : samples) {
            if (sample.isSomaticTumorSample()) {
                tumorSamples.add(sample);
            }
        }
        return tumorSamples;
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

    private static boolean checkInputs(@NotNull CommandLine cmd) {
        String runsDirectory = cmd.getOptionValue(RUNS_DIRECTORY);

        boolean allParamsPresent = !Utils.anyNull(runsDirectory,
                cmd.getOptionValue(CPCT_ECRF_FILE),
                cmd.getOptionValue(CPCT_FORM_STATUS_CSV),
                cmd.getOptionValue(DRUP_ECRF_FILE),
                cmd.getOptionValue(LIMS_DIRECTORY),
                cmd.getOptionValue(TREATMENT_MAPPING_CSV),
                cmd.getOptionValue(BIOPSY_MAPPING_CSV),
                cmd.getOptionValue(TUMOR_LOCATION_MAPPING_TSV),
                cmd.getOptionValue(CURATED_PRIMARY_TUMOR_TSV),
                cmd.getOptionValue(PATIENT_TUMOR_CURATION_STATUS_TSV),
                cmd.getOptionValue(DOID_JSON),
                cmd.getOptionValue(REPORTING_DB_TSV));

        if (cmd.hasOption(DO_LOAD_CLINICAL_DATA)) {
            allParamsPresent = allParamsPresent && DatabaseAccess.hasDatabaseConfig(cmd);
        }

        if (cmd.hasOption(DO_PROCESS_WIDE_CLINICAL_DATA)) {
            allParamsPresent = allParamsPresent && !Utils.anyNull(cmd.getOptionValue(WIDE_AVL_TREATMENT_CSV),
                    cmd.getOptionValue(WIDE_PRE_AVL_TREATMENT_CSV),
                    cmd.getOptionValue(WIDE_BIOPSY_CSV),
                    cmd.getOptionValue(WIDE_RESPONSE_CSV),
                    cmd.getOptionValue(WIDE_FIVE_DAYS_CSV));
        }

        boolean validRunDirectories = true;
        if (allParamsPresent) {
            File runDirectoryDb = new File(runsDirectory);

            if (!runDirectoryDb.exists() || !runDirectoryDb.isDirectory()) {
                validRunDirectories = false;
                LOGGER.warn("HMF database run directory '{}' does not exist or is not a directory", runDirectoryDb);
            }
        }

        return validRunDirectories && allParamsPresent;
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(RUNS_DIRECTORY,
                true,
                "Path towards the folder containing patient runs that are considered part of HMF database.");

        options.addOption(REPORTING_DB_TSV, true, "Path towards the reporting db tsv file.");

        options.addOption(CPCT_ECRF_FILE, true, "Path towards the CPCT ecrf file.");
        options.addOption(CPCT_FORM_STATUS_CSV, true, "Path towards the CPCT form status csv file.");
        options.addOption(DRUP_ECRF_FILE, true, "Path towards the DRUP ecrf file.");
        options.addOption(DO_LOAD_RAW_ECRF, false, "If set, writes raw ecrf data to database.");

        options.addOption(DO_LOAD_CLINICAL_DATA, false, "If set, clinical data will be loaded into the database.");
        options.addOption(CURATED_PRIMARY_TUMOR_TSV, true, "Path towards to the curated primary tumor TSV.");
        options.addOption(PATIENT_TUMOR_CURATION_STATUS_TSV, true, "Path where patient tumor curation status will be written to");

        options.addOption(DO_PROCESS_WIDE_CLINICAL_DATA,
                false,
                "if set, creates clinical timeline for WIDE patients and persists to database.");
        options.addOption(WIDE_AVL_TREATMENT_CSV, true, "Path towards the WIDE avl treatment csv");
        options.addOption(WIDE_PRE_AVL_TREATMENT_CSV, true, "Path towards the WIDE pre avl treatment csv.");
        options.addOption(WIDE_BIOPSY_CSV, true, "Path towards the WIDE biopsy csv.");
        options.addOption(WIDE_RESPONSE_CSV, true, "Path towards the WIDE response csv.");
        options.addOption(WIDE_FIVE_DAYS_CSV, true, "Path towards the WIDE five days csv.");

        options.addOption(LIMS_DIRECTORY, true, "Path towards the LIMS directory.");

        options.addOption(DOID_JSON, true, "Path towards to the json file of the doid ID of primary tumors.");
        options.addOption(TUMOR_LOCATION_MAPPING_TSV, true, "Path towards to the tumor location mapping TSV.");
        options.addOption(TREATMENT_MAPPING_CSV, true, "Path towards to the treatment mapping CSV.");
        options.addOption(BIOPSY_MAPPING_CSV, true, "Path towards to the biopsy mapping CSV.");

        options.addOption(PIPELINE_VERSION, true, "Path towards the pipeline version");

        addDatabaseCmdLineArgs(options);
        return options;
    }
}
