package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
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
    private static final String DB_USER = "dbUser";
    private static final String DB_PASS = "dbPass";
    private static final String DB_URL = "dbUrl";

    public static void main(String[] args)
            throws ParseException, IOException, InterruptedException, java.text.ParseException, XMLStreamException,
            SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String runFolderPath = cmd.getOptionValue(RUNS_DIR);
        final String ecrfFilePath = cmd.getOptionValue(ECRF_FILE);
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        if (runFolderPath == null || ecrfFilePath == null || userName == null || password == null
                || databaseUrl == null) {
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
                final CpctPatientReader cpctPatientReader = new CpctPatientReader(model);
                final DatabaseWriter dbWriter = new DatabaseWriter(userName, password, jdbcUrl);
                dbWriter.clearTables();
                for (final EcrfPatient patient : patients) {
                    final Patient cpctPatient = cpctPatientReader.read(patient);
                    dbWriter.writePatient(cpctPatient);
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
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}