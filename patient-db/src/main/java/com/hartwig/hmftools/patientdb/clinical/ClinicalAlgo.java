package com.hartwig.hmftools.patientdb.clinical;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfig;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfigFactory;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.EcrfModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.lims.Lims;
import com.hartwig.hmftools.patientdb.clinical.readers.LimsPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.ColoPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.CorePatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.EcrfPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctPatientReader;
import com.hartwig.hmftools.patientdb.clinical.readers.cpct.CpctUtil;
import com.hartwig.hmftools.patientdb.clinical.readers.drup.DrupPatientReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import static java.lang.String.join;

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
    @NotNull
    private final List<DoidNode> doidNodes;

    ClinicalAlgo(@NotNull final EcrfModels ecrfModels, @NotNull final PrimaryTumorCurator primaryTumorCurator,
            @NotNull final BiopsySiteCurator biopsySiteCurator, @NotNull final TreatmentCurator treatmentCurator, List<DoidNode> doidNodes) {
        this.ecrfModels = ecrfModels;
        this.primaryTumorCurator = primaryTumorCurator;
        this.biopsySiteCurator = biopsySiteCurator;
        this.treatmentCurator = treatmentCurator;
        this.doidNodes = doidNodes;
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

        LOGGER.info("Interpreting and curating data for LIMS patients");
        List<Patient> limsPatients = readLimsPatients(samplesPerPatient, consentConfigMap, doidNodes);
        LOGGER.info(" Finished curation of {} LIMS patients", limsPatients.size());

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

        List<Patient> mergedPatients = Lists.newArrayList();
        mergedPatients.addAll(readColoPatients(lims));
        mergedPatients.addAll(limsPatients);
        mergedPatients.addAll(cpctPatients);
        mergedPatients.addAll(drupPatients);

        return mergedPatients;
    }

    @NotNull
    private static List<Patient> readEcrfPatients(@NotNull EcrfPatientReader reader, @NotNull Iterable<EcrfPatient> ecrfPatients,
            @NotNull Map<String, List<SampleData>> samplesPerPatient, @NotNull Map<String, ConsentConfig> consentConfigMap)
            throws IOException {
        List<Patient> patients = Lists.newArrayList();
        for (EcrfPatient ecrfPatient : ecrfPatients) {
            List<SampleData> sequencedSamples = sequencedSamplesOnly(samplesPerPatient.get(ecrfPatient.patientId()));
            String cohortId = !sequencedSamples.isEmpty() ? sequencedSamples.get(0).cohortId() : Strings.EMPTY;

            patients.add(reader.read(ecrfPatient, sequencedSamples, consentConfigMap, cohortId));
        }
        return patients;
    }

    @NotNull
    private List<Patient> readLimsPatients(@NotNull Map<String, List<SampleData>> samplesPerPatient,
            Map<String, ConsentConfig> consentConfigMap, List<DoidNode> doidNodes) {
        List<Patient> patients = Lists.newArrayList();
        CorePatientReader corePatientReader = new CorePatientReader(primaryTumorCurator);
        LimsPatientReader limsPatientReader = new LimsPatientReader(doidNodes);

        for (Map.Entry<String, List<SampleData>> entry : samplesPerPatient.entrySet()) {
            List<SampleData> tumorSamples = tumorSamplesOnly(entry.getValue());
            if (!tumorSamples.isEmpty()) {
                String patientId = entry.getKey();

                // We assume every sample for a single patient has the same primary tumor.
                SampleData chosenSample = tumorSamples.get(0);
                String logInfo = join("|", chosenSample.sampleBarcode(), chosenSample.sampleId(),
                        chosenSample.cohortId());

                if (!hasAllRequiredFields(chosenSample)) {
                    LOGGER.warn("Not processing sample from LIMS/LAMA due to required fields missing [{}]", logInfo);
                } else if (hasLamaBasedCuration(chosenSample)) {
                    LOGGER.info("Processing sample from LIMS/LAMA sample with curation from LAMA [{}]", logInfo);
                    patients.add(limsPatientReader.read(patientId, chosenSample, sequencedSamplesOnly(tumorSamples)));
                } else if (isCuratedFromEcrfSource(chosenSample)) {
                    LOGGER.info("Not processing sample from LIMS/LAMA because has ECRF cohort [{}]", logInfo);
                } else if (hasPrimaryTumorStringForCuration(chosenSample)) {
                    logInfo = join("|", logInfo, "ptum:" + chosenSample.limsPrimaryTumor());
                    LOGGER.info("Processing sample from LIMS/LAMA without curation but with ptum string [{}]", logInfo);
                    patients.add(corePatientReader.read(patientId,
                            chosenSample.limsPrimaryTumor(),
                            sequencedSamplesOnly(tumorSamples),
                            consentConfigMap,
                            chosenSample.cohortId()));
                } else {
                    LOGGER.warn("Not processing from LIMS/LAMA because no tumor information present [{}]", logInfo);
                }
            }
        }
        return patients;
    }

    private boolean isCuratedFromEcrfSource(SampleData sample) {
        String cohortStart = sample.cohortId().substring(0, 4);
        return List.of("CPCT", "DRUP").contains(cohortStart);
    }

    private boolean hasAllRequiredFields(SampleData sample) {
        return !Objects.requireNonNull(sample.cohortId()).isEmpty();
    }

    private boolean hasLamaBasedCuration(SampleData sample) {
        return sample.limsTumorLocation() != null && !Objects.requireNonNull(sample.limsTumorLocation()).isEmpty();
    }

    private boolean hasPrimaryTumorStringForCuration(SampleData sample) {
        return sample.limsPrimaryTumor() != null && !Objects.requireNonNull(sample.limsPrimaryTumor()).isEmpty();
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
