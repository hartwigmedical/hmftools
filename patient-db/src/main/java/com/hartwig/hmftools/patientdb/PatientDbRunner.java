package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PatientDbRunner {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final String RUNS_DIR = "runsDir";
    private static final String ECRF_FILE = "ecrf";

    public static void main(String[] args)
            throws ParseException, IOException, InterruptedException, java.text.ParseException, XMLStreamException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String runFolderPath = cmd.getOptionValue(RUNS_DIR);
        final String ecrfFilePath = cmd.getOptionValue(ECRF_FILE);

        if (runFolderPath == null) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-Db", options);
        } else {
            final File dir = new File(runFolderPath);
            if (dir.isDirectory()) {
                final List<CpctRunData> data = RunsFolderProcessor.getPatientRunsData(dir);
                LOGGER.info("Listing data for " + data.size() + " patients.");
                LOGGER.info(data.toString());
                LOGGER.info("Loading ecrf model...");
                final CpctEcrfModel model = CpctEcrfModel.loadFromXML(ecrfFilePath);
                final List<String> cpctPatientIds = data.stream().map(CpctRunData::patientId).filter(
                        id -> id.startsWith("CPCT")).collect(Collectors.toList());
                final Iterable<EcrfPatient> patients = model.findPatientsById(cpctPatientIds);
                LOGGER.info("Reading CPCT patient data...");
                final CpctPatientDataReader cpctPatientDataReader = new CpctPatientDataReader(model);
                final CpctTumorDataReader cpctTumorDataReader = new CpctTumorDataReader();
                final CpctSystemicTherapyReader cpctSystemicTherapyReader = new CpctSystemicTherapyReader();
                final CpctRadioTherapyReader cpctRadioTherapyReader = new CpctRadioTherapyReader();
                final CpctTreatmentDataReader cpctTreatmentDataReader = new CpctTreatmentDataReader();
                for (final EcrfPatient patient : patients) {
                    final PatientData patientData = cpctPatientDataReader.read(patient);
                    final Optional<TumorData> tumorDataOpt = cpctTumorDataReader.read(patient);
                    final Optional<List<SystemicTherapyData>> systemicTherapiesOpt = cpctSystemicTherapyReader.read(
                            patient);
                    final Optional<List<RadioTherapyData>> radioTherapiesOpt = cpctRadioTherapyReader.read(patient);
                    final Optional<TreatmentData> treatmentDataOpt = cpctTreatmentDataReader.read(patient);
                    LOGGER.info(patientData.toString());
                    tumorDataOpt.ifPresent(tumorData -> LOGGER.info(tumorData.toString()));
                    systemicTherapiesOpt.ifPresent(systemicTherapies -> LOGGER.info(systemicTherapies.toString()));
                    radioTherapiesOpt.ifPresent(radioTherapies -> LOGGER.info(radioTherapies.toString()));
                    treatmentDataOpt.ifPresent(treatmentData -> LOGGER.info(treatmentData.toString()));
                }
            } else {
                if (!dir.exists()) {
                    LOGGER.warn("dir " + dir + " does not exist.");
                }
                HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("Patient-Db", options);
            }
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(RUNS_DIR, true, "Path towards the folder containing patient runs.");
        options.addOption(ECRF_FILE, true, "Path towards the cpct ecrf file.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}