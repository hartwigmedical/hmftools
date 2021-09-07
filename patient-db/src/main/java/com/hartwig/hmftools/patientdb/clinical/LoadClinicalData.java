package com.hartwig.hmftools.patientdb.clinical;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_URL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.EcrfModel;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;
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
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LoadClinicalData {

    private static final Logger LOGGER = LogManager.getLogger(LoadClinicalData.class);
    private static final String VERSION = LoadClinicalData.class.getPackage().getImplementationVersion();

    public static void main(@NotNull String[] args) throws IOException, XMLStreamException, SQLException, ParseException {
        LOGGER.info("Running Clinical Patient DB v{}", VERSION);
        Options options = ClinicalAlgoConfig.createOptions();

        ClinicalAlgoConfig config = null;
        try {
            config = ClinicalAlgoConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("Clinical-Patient-DB", options);
            System.exit(1);
        }

        LOGGER.info("Creating lims from {}", config.limsDirectory());
        Lims lims = LimsFactory.fromLimsDirectory(config.limsDirectory());

        LOGGER.info("Creating clinical algorithm");
        ClinicalAlgo algo = ClinicalAlgoBuilder.fromConfig(config);

        Map<String, List<SampleData>> samplesPerPatient = new SampleLoader(lims).loadSamplesPerPatient(config);

        LOGGER.info("Running clinical patient db algo on {} patients", samplesPerPatient.size());
        List<Patient> patients = algo.interpret(samplesPerPatient, lims);

        LOGGER.info("Writing primary tumor data for {} patients", samplesPerPatient.size());
        PrimaryTumorDataWriter.write(config, samplesPerPatient, patients);

        if (config.doLoadClinicalData()) {
            DatabaseAccess.addDatabaseCmdLineArgs(options);
            CommandLine cmd = new DefaultParser().parse(options, args);

            LOGGER.info("Writing clinical data to database '{}'", cmd.getOptionValue(DB_URL));
            DatabaseAccess dbWriter = databaseAccess(cmd);

            if (config.doLoadRawEcrf()) {
                writeRawEcrf(dbWriter, algo.ecrfModels(), samplesPerPatient);
            }

            writeClinicalData(dbWriter, lims, samplesPerPatient, patients);

            dbWriter.writeValidationFindings(CurationValidator.validatePrimaryTumorCurator(algo.primaryTumorCurator()));
            dbWriter.writeValidationFindings(CurationValidator.validateTreatmentCurator(algo.treatmentCurator()));
        }

        LOGGER.info("Complete");
    }

    private static void writeRawEcrf(@NotNull DatabaseAccess dbWriter, @NotNull EcrfModels ecrfModels,
            @NotNull Map<String, List<SampleData>> samplesPerPatient) {
        Set<String> sequencedPatients = sequencedPatients(samplesPerPatient);

        EcrfModel cpctEcrfModel = ecrfModels.cpctModel();
        LOGGER.info(" Writing raw cpct ecrf data for {} patients", cpctEcrfModel.patientCount());
        dbWriter.clearCpctEcrf();
        dbWriter.writeCpctEcrf(cpctEcrfModel, sequencedPatients);
        LOGGER.info("  Finished writing raw cpct ecrf data for {} patients", cpctEcrfModel.patientCount());

        EcrfModel drupEcrfModel = ecrfModels.drupModel();
        LOGGER.info(" Writing raw drup ecrf data for {} patients", drupEcrfModel.patientCount());
        dbWriter.clearDrupEcrf();
        dbWriter.writeDrupEcrf(drupEcrfModel, sequencedPatients);
        LOGGER.info("  Finished writing raw drup ecrf data for {} patients", drupEcrfModel.patientCount());
    }

    private static void writeClinicalData(@NotNull DatabaseAccess dbAccess, @NotNull Lims lims,
            @NotNull Map<String, List<SampleData>> samplesPerPatient, @NotNull List<Patient> patients) {
        LOGGER.info(" Clearing interpreted clinical tables in database");
        dbAccess.clearClinicalTables();

        Set<String> sequencedPatients = sequencedPatients(samplesPerPatient);

        int missingPatients = 0;
        int missingSamples = 0;
        LOGGER.info(" Writing clinical data for {} sequenced patients", sequencedPatients.size());
        for (String patientId : sequencedPatients) {
            Patient patient = findByPatientId(patients, patientId);
            if (patient == null) {
                LOGGER.warn(" No clinical data found for patient {}", patientId);
                missingPatients++;
                List<SampleData> sequencedSamples = sequencedSamples(samplesPerPatient.get(patientId));
                missingSamples += sequencedSamples.size();
                dbAccess.writeSampleClinicalData(patientId, lims.isBlacklisted(patientId), sequencedSamples);
            } else if (patient.sequencedBiopsies().isEmpty()) {
                LOGGER.warn(" No sequenced biopsies found for sequenced patient: {}! Skipping writing to db", patientId);
            } else {
                dbAccess.writeFullClinicalData(patient, lims.isBlacklisted(patientId));
                List<ValidationFinding> findings = PatientValidator.validatePatient(patient);

                dbAccess.writeValidationFindings(findings);
                dbAccess.writeValidationFindings(patient.matchFindings());
            }
        }

        if (missingPatients > 0) {
            LOGGER.warn(" Could not load {} patients ({} samples)!", missingPatients, missingSamples);
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

    @NotNull
    private static Set<String> sequencedPatients(@NotNull Map<String, List<SampleData>> samplesPerPatient) {
        Set<String> sequencedPatients = Sets.newHashSet();
        for (Map.Entry<String, List<SampleData>> entry : samplesPerPatient.entrySet()) {
            for (SampleData sample : entry.getValue()) {
                if (sample.sequenced()) {
                    sequencedPatients.add(entry.getKey());
                }
            }
        }
        return sequencedPatients;
    }

    @NotNull
    private static List<SampleData> sequencedSamples(@Nullable Iterable<SampleData> samples) {
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
