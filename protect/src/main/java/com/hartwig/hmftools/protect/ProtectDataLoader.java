package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_PASS;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_URL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_USER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.protect.actionability.EvidenceItem;
import com.hartwig.hmftools.protect.common.ActionabilityFile;


import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProtectDataLoader {
    private static final Logger LOGGER = LogManager.getLogger(ProtectDataLoader.class);

    private static final String SAMPLE = "sample";

    private static final String ACTIONABILITY_TSV = "knowledgebase_dir";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        final String sampleId = cmd.getOptionValue(SAMPLE);

        final String actionabilityTsv = cmd.getOptionValue(ACTIONABILITY_TSV);

        if (anyNull(sampleId,
                actionabilityTsv,
                cmd.getOptionValue(DB_USER),
                cmd.getOptionValue(DB_PASS),
                cmd.getOptionValue(DB_URL))) {
            printUsageAndExit(options);
        }

        LOGGER.info("Connecting with database");
        DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Reading actionability for sample {}", sampleId);
        List<EvidenceItem> combinedEvidence = ActionabilityFile.read(actionabilityTsv);

        LOGGER.info("Writing evidence items into db");
        final ClinicalEvidenceDAOProtect clinicalEvidenceDAOProtect = new ClinicalEvidenceDAOProtect(dbAccess.context());
        clinicalEvidenceDAOProtect.writeClinicalEvidence(sampleId, combinedEvidence);
        LOGGER.info("Finished");
    }

    public static boolean anyNull(@NotNull final Object... arguments) {
        for (final Object object : arguments) {
            if (object == null) {
                return true;
            }
        }
        return false;
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("protect - evidence data loader", options);
        System.exit(1);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(SAMPLE, true, "Tumor sample of run");

        options.addOption(ACTIONABILITY_TSV, true, "Path towards the TSV file of the actionability TSV");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}