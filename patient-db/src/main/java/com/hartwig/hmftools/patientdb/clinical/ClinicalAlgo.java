package com.hartwig.hmftools.patientdb.clinical;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfig;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfigFactory;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.EcrfModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.readers.ColoPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.CorePatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.EcrfPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.WidePatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctUtil;
import com.hartwig.hmftools.patientdb.clinical.readers.drup.DrupPatientReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ClinicalAlgo {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalAlgo.class);

    @NotNull
    private final EcrfModels ecrfModels;
    @NotNull
    private final PrimaryTumorCurator primaryTumorCurator;
    @NotNull
    private final BiopsySiteCurator biopsySiteCurator;
    @NotNull
    private final TreatmentCurator treatmentCurator;

    ClinicalAlgo(@NotNull final EcrfModels ecrfModels, @NotNull final PrimaryTumorCurator primaryTumorCurator,
            @NotNull final BiopsySiteCurator biopsySiteCurator, @NotNull final TreatmentCurator treatmentCurator) {
        this.ecrfModels = ecrfModels;
        this.primaryTumorCurator = primaryTumorCurator;
        this.biopsySiteCurator = biopsySiteCurator;
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    public EcrfModels ecrfModels() {
        return ecrfModels;
    }

    @NotNull
    public PrimaryTumorCurator primaryTumorCurator() {
        return primaryTumorCurator;
    }

    @NotNull
    public TreatmentCurator treatmentCurator() {
        return treatmentCurator;
    }

    @NotNull
    public List<Patient> interpret(@NotNull Map<String, List<SampleData>> samplesPerPatient, @NotNull Lims lims,
            @NotNull String consentConfigTsv) throws IOException {
        EcrfModel cpctEcrfModel = ecrfModels.cpctModel();

        LOGGER.info("Reading the informed consent configuration file");
        Map<String, ConsentConfig> consentConfigMap = ConsentConfigFactory.read(consentConfigTsv);

        LOGGER.info("Interpreting and curating data for {} CPCT patients", cpctEcrfModel.patientCount());
        EcrfPatientReader cpctPatientReader =
                new CpctPatientReader(primaryTumorCurator, CpctUtil.extractHospitalMap(cpctEcrfModel), biopsySiteCurator, treatmentCurator);

        List<Patient> cpctPatients = readEcrfPatients(cpctPatientReader, cpctEcrfModel.patients(), samplesPerPatient, consentConfigMap);
        LOGGER.info(" Finished curation of {} CPCT patients", cpctPatients.size());

        EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info("Interpreting and curating data for {} DRUP patients", drupEcrfModel.patientCount());
        EcrfPatientReader drupPatientReader = new DrupPatientReader(primaryTumorCurator, biopsySiteCurator);

        List<Patient> drupPatients = readEcrfPatients(drupPatientReader, drupEcrfModel.patients(), samplesPerPatient, consentConfigMap);
        LOGGER.info(" Finished curation of {} DRUP patients", drupPatients.size());

        LOGGER.info("Interpreting and curating data for WIDE patients");
        List<Patient> widePatients = readWidePatients(samplesPerPatient, consentConfigMap);
        LOGGER.info(" Finished curation of {} WIDE patients", widePatients.size());

        LOGGER.info("Interpreting and curating data for CORE patients");
        List<Patient> corePatients = readCoreAndActinPatients(samplesPerPatient, "CORE", consentConfigMap);
        LOGGER.info(" Finished curation of {} CORE patients", corePatients.size());

        LOGGER.info("Interpreting and curating data for ACTIN patients");
        List<Patient> actinPatients = readCoreAndActinPatients(samplesPerPatient, "ACTIN", consentConfigMap);
        LOGGER.info(" Finished curation of {} ACTIN patients", actinPatients.size());

        List<Patient> mergedPatients = Lists.newArrayList();
        mergedPatients.addAll(cpctPatients);
        mergedPatients.addAll(drupPatients);
        mergedPatients.addAll(widePatients);
        mergedPatients.addAll(corePatients);
        mergedPatients.addAll(actinPatients);
        mergedPatients.addAll(readColoPatients(lims));
        return mergedPatients;
    }

    @NotNull
    private static List<Patient> readEcrfPatients(@NotNull EcrfPatientReader reader, @NotNull Iterable<EcrfPatient> ecrfPatients,
            @NotNull Map<String, List<SampleData>> samplesPerPatient, @NotNull Map<String, ConsentConfig> consentConfigMap)
            throws IOException {
        List<Patient> patients = Lists.newArrayList();

        for (EcrfPatient ecrfPatient : ecrfPatients) {
            List<SampleData> sequencedSamples = sequencedSamplesOnly(samplesPerPatient.get(ecrfPatient.patientId()));
            patients.add(reader.read(ecrfPatient, sequencedSamples, consentConfigMap, sequencedSamples.get(0).cohortId()));
        }
        return patients;
    }

    @NotNull
    private List<Patient> readWidePatients(@NotNull Map<String, List<SampleData>> samplesPerPatient,
            @NotNull Map<String, ConsentConfig> consentConfigMap) {
        List<Patient> patients = Lists.newArrayList();

        WidePatientReader widePatientReader = new WidePatientReader(ecrfModels.wideModel(), primaryTumorCurator, treatmentCurator);
        for (Map.Entry<String, List<SampleData>> entry : samplesPerPatient.entrySet()) {
            List<SampleData> tumorSamples = tumorSamplesOnly(entry.getValue());
            if (!tumorSamples.isEmpty() && tumorSamples.get(0).cohortId().equals("WIDE")) {
                String patientId = entry.getKey();
                // We assume every sample for a single patient has the same primary tumor.
                String primaryTumor = tumorSamples.get(0).limsPrimaryTumor();
                patients.add(widePatientReader.read(patientId,
                        primaryTumor,
                        sequencedSamplesOnly(tumorSamples),
                        consentConfigMap,
                        tumorSamples.get(0).cohortId()));
            }
        }
        return patients;
    }

    @NotNull
    private List<Patient> readCoreAndActinPatients(@NotNull Map<String, List<SampleData>> samplesPerPatient, @NotNull String cohort,
            Map<String, ConsentConfig> consentConfigMap) {
        List<Patient> patients = Lists.newArrayList();
        CorePatientReader corePatientReader = new CorePatientReader(primaryTumorCurator);

        for (Map.Entry<String, List<SampleData>> entry : samplesPerPatient.entrySet()) {
            List<SampleData> tumorSamples = tumorSamplesOnly(entry.getValue());
            if (!tumorSamples.isEmpty() && tumorSamples.get(0).cohortId().contains(cohort)) {
                String patientId = entry.getKey();
                // We assume every sample for a single patient has the same primary tumor.
                String primaryTumor = tumorSamples.get(0).limsPrimaryTumor();
                patients.add(corePatientReader.read(patientId,
                        primaryTumor,
                        sequencedSamplesOnly(tumorSamples),
                        consentConfigMap,
                        tumorSamples.get(0).cohortId()));
            }
        }

        return patients;
    }

    @NotNull
    private List<Patient> readColoPatients(@NotNull Lims lims) {
        List<Patient> patients = Lists.newArrayList();

        String tumorLocation = lims.primaryTumor("COLO829V003TVAL");
        String patientId = lims.patientId("COLO829V003TVAL");

        ColoPatientReader coloPatientReader = new ColoPatientReader(primaryTumorCurator);
        LOGGER.info("Creating patient representation for COLO829");
        patients.add(coloPatientReader.read(patientId, tumorLocation));

        return patients;
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
    private static List<SampleData> sequencedSamplesOnly(@Nullable Iterable<SampleData> samples) {
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
}
