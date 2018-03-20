package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.formstatus.ImmutableFormStatusModel;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.readers.RunsFolderReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class LoadDrupEcrfData {
    private static final Logger LOGGER = LogManager.getLogger(LoadDrupEcrfData.class);

    private static final String RUNS_DIR = "runs_dir";
    private static final String ECRF_FILE = "ecrf";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, XMLStreamException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);

        final String ecrfFile = cmd.getOptionValue(ECRF_FILE);
        final String runsFolderPath = cmd.getOptionValue(RUNS_DIR);

        if (Utils.anyNull(userName, password, databaseUrl, ecrfFile, runsFolderPath)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load DRUP ecrf", options);
        } else {
            final File runsDirectory = new File(runsFolderPath);
            if (runsDirectory.isDirectory()) {
                final String jdbcUrl = "jdbc:" + databaseUrl;
                final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);
                dbWriter.clearDrupEcrf();
                LOGGER.info("Importing DRUP ecrf data from: {}", ecrfFile);
                final CpctEcrfModel model = CpctEcrfModel.loadFromXML(ecrfFile, new ImmutableFormStatusModel(Maps.newHashMap()));
                final List<RunContext> runContexts = RunsFolderReader.getRunContexts(runsDirectory);
                final Set<String> sequencedPatients = Utils.sequencedPatientIds(runContexts);
                LOGGER.info("Writing raw ecrf data for " + model.patientCount() + " patients.");
                dbWriter.writeDrupEcrf(model, sequencedPatients);
                LOGGER.info("Done writing raw ecrf data for " + model.patientCount() + " patients!");
            } else {
                if (!runsDirectory.exists()) {
                    LOGGER.warn("dir " + runsDirectory + " does not exist.");
                }
                final HelpFormatter formatter = new HelpFormatter();
                formatter.printHelp("patient-db - load DRUP ecrf", options);
            }
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(RUNS_DIR, true, "Path towards the folder containing patient runs.");
        options.addOption(ECRF_FILE, true, "Path towards the DRUP ecrf file.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
