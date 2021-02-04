package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_URL;
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
import com.hartwig.hmftools.patientdb.clinical.DumpPrimaryTumorData;
import com.hartwig.hmftools.patientdb.clinical.EcrfModels;
import com.hartwig.hmftools.patientdb.clinical.ImmutableEcrfModels;
import com.hartwig.hmftools.patientdb.clinical.LoadClinicalDataConfig;
import com.hartwig.hmftools.patientdb.clinical.context.RunContext;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.readers.ColoPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.CorePatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.EcrfPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.LimsSampleReader;
import com.hartwig.hmftools.patientdb.clinical.readers.RunsFolderReader;
import com.hartwig.hmftools.patientdb.clinical.readers.WidePatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctUtil;
import com.hartwig.hmftools.patientdb.clinical.readers.drup.DrupPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.ImmutableWideEcrfModel;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideAvlTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideBiopsyData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideEcrfFileReader;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideEcrfModel;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideFiveDays;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WidePreAvlTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.readers.wide.WideResponseData;
import com.hartwig.hmftools.patientdb.clinical.validators.CurationValidator;
import com.hartwig.hmftools.patientdb.clinical.validators.PatientValidator;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LoadClinicalData {

    private static final Logger LOGGER = LogManager.getLogger(LoadClinicalData.class);
    private static final String VERSION = LoadClinicalData.class.getPackage().getImplementationVersion();

    public static void main(@NotNull String[] args) throws IOException, XMLStreamException, SQLException, ParseException {
        LOGGER.info("Running Clinical Patient DB v{}", VERSION);
        Options options = LoadClinicalDataConfig.createOptions();

        LoadClinicalDataConfig config = null;
        try {
            config = LoadClinicalDataConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("Clinical-Patient-DB", options);
            System.exit(1);
        }

        List<DoidNode> doidNodes = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJson()).nodes();
        PrimaryTumorCurator primaryTumorCurator = new PrimaryTumorCurator(config.tumorLocationMappingTsv(), doidNodes);
        BiopsySiteCurator biopsySiteCurator = new BiopsySiteCurator(config.biopsyMappingCsv());
        TreatmentCurator treatmentCurator = new TreatmentCurator(config.treatmentMappingCsv());

        LOGGER.info("Loading sequence runs from {}", config.runsDirectory());
        List<RunContext> runContexts = loadRunContexts(config.runsDirectory(), config.pipelineVersionFile());
        Map<String, List<String>> sequencedSamplesPerPatient = extractSequencedSamplesFromRunContexts(runContexts);
        Map<String, String> sampleToSetNameMap = extractSampleToSetNameMap(runContexts);
        Set<String> sequencedPatientIds = sequencedSamplesPerPatient.keySet();

        LOGGER.info(" Loaded sequence runs for {} patient IDs ({} samples)",
                sequencedPatientIds.size(),
                toUniqueSampleIds(sequencedSamplesPerPatient).size());

        LOGGER.info("Loading sample data from LIMS in {}", config.limsDirectory());
        Lims lims = LimsFactory.fromLimsDirectory(config.limsDirectory());
        Map<String, List<SampleData>> sampleDataPerPatient =
                extractAllSamplesFromLims(lims, sampleToSetNameMap, sequencedSamplesPerPatient);
        LOGGER.info(" Loaded samples for {} patient IDs ({} samples)",
                sampleDataPerPatient.keySet().size(),
                countValues(sampleDataPerPatient));

        EcrfModels ecrfModels = loadEcrfModels(config);

        List<Patient> patients = interpret(sampleDataPerPatient, ecrfModels, primaryTumorCurator, biopsySiteCurator, treatmentCurator);

        LOGGER.info("Check for missing curation tumor location when info is known");
        Map<String, PatientTumorCurationStatus> patientTumorCurationStatusMap =
                generatePatientTumorCurationStatusMap(sampleDataPerPatient, patients, config.reportingDbTsv());

        LOGGER.info("Writing patient tumor curation status");
        DumpPrimaryTumorData.writePatientTumorCurationStatesToTSV(config.patientTumorCurationStatusTsv(), patientTumorCurationStatusMap);

        LOGGER.info("Writing curated primary tumors");
        DumpPrimaryTumorData.writeCuratedPrimaryTumorsToTSV(config.curatedPrimaryTumorTsv(), patients);

        if (config.doLoadClinicalData()) {
            DatabaseAccess.addDatabaseCmdLineArgs(options);
            CommandLine cmd = new DefaultParser().parse(options, args);

            LOGGER.info("Connecting to database {}", cmd.getOptionValue(DB_URL));
            DatabaseAccess dbWriter = databaseAccess(cmd);

            if (config.doLoadRawEcrf()) {
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
                if (!patientId.startsWith("CORE") || !patientId.startsWith("WIDE")) {
                    patientTumorCurationStatusMap.put(patientId, PatientTumorCurationStatus.NOT_RESOLVED);
                }
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
            String patientId = extractPatientIdentifier(runContext.setName());
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
    private static EcrfModels loadEcrfModels(@NotNull LoadClinicalDataConfig config) throws IOException, XMLStreamException {
        EcrfModel cpctEcrfModel = buildCpctEcrfModel(config);
        EcrfModel drupEcrfModel = buildDrupEcrfModel(config);
        WideEcrfModel wideEcrfModel = buildWideEcrfModel(config);

        return ImmutableEcrfModels.builder().cpctModel(cpctEcrfModel).drupModel(drupEcrfModel).wideModel(wideEcrfModel).build();
    }

    @NotNull
    private static EcrfModel buildCpctEcrfModel(@NotNull LoadClinicalDataConfig config) throws IOException, XMLStreamException {
        LOGGER.info("Loading CPCT eCRF from {}", config.cpctEcrfFile());
        FormStatusModel cpctFormStatusModel = FormStatusReader.buildModelFromCsv(config.cpctFormStatusCsv());
        EcrfModel cpctEcrfModel = EcrfModel.loadFromXMLWithFormStates(config.cpctEcrfFile(), cpctFormStatusModel);
        LOGGER.info(" Finished loading CPCT eCRF. Read {} patients", cpctEcrfModel.patientCount());

        return cpctEcrfModel;
    }

    @NotNull
    private static EcrfModel buildDrupEcrfModel(@NotNull LoadClinicalDataConfig config) throws FileNotFoundException, XMLStreamException {
        LOGGER.info("Loading DRUP eCRF from {}", config.drupEcrfFile());
        EcrfModel drupEcrfModel = EcrfModel.loadFromXMLNoFormStates(config.drupEcrfFile());
        LOGGER.info(" Finished loading DRUP eCRF. Read {} patients", drupEcrfModel.patientCount());

        return drupEcrfModel;
    }

    @NotNull
    private static WideEcrfModel buildWideEcrfModel(@NotNull LoadClinicalDataConfig config) throws IOException {
        WideEcrfModel wideEcrfModel;

        if (config.doProcessWideClinicalData()) {
            LOGGER.info("Loading WIDE eCRF");

            String preAvlTreatmentCsv = config.widePreAvlTreatmentCsv();
            List<WidePreAvlTreatmentData> preAvlTreatments = WideEcrfFileReader.readPreAvlTreatments(preAvlTreatmentCsv);
            LOGGER.info(" Loaded {} WIDE pre-AVL-treatments from {}", preAvlTreatments.size(), preAvlTreatmentCsv);

            String biopsyCsv = config.wideBiopsyCsv();
            List<WideBiopsyData> biopsies = WideEcrfFileReader.readBiopsies(biopsyCsv);
            LOGGER.info(" Loaded {} WIDE biopsies from {}", biopsies.size(), biopsyCsv);

            String avlTreatmentCsv = config.wideAvlTreatmentCsv();
            List<WideAvlTreatmentData> avlTreatments = WideEcrfFileReader.readAvlTreatments(avlTreatmentCsv);
            LOGGER.info(" Loaded {} WIDE AVL treatments from {}", avlTreatments.size(), avlTreatmentCsv);

            String wideResponseCsv = config.wideResponseCsv();
            List<WideResponseData> responses = WideEcrfFileReader.readResponses(wideResponseCsv);
            LOGGER.info(" Loaded {} WIDE responses from {}", responses.size(), wideResponseCsv);

            String fiveDaysCsv = config.wideFiveDaysCsv();
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

    @NotNull
    private static String extractPatientIdentifier(@NotNull String setName) {
        String[] names = setName.split("_");
        if (names.length < 5) {
            LOGGER.error("Run name {} had less than 5 parts after splitting on _", setName);
            return Strings.EMPTY;
        }
        return names[4];
    }
}
