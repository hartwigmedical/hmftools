package com.hartwig.hmftools.iclusion.dao;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadIclusionDatabase {

    private static final Logger LOGGER = LogManager.getLogger(LoadIclusionDatabase.class);

    private static final String ICLUSION_DB_TSV = "iclusion_db_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String iclusionDBTsv = cmd.getOptionValue(ICLUSION_DB_TSV);

        if (Utils.anyNull(iclusionDBTsv) || !new File(iclusionDBTsv).exists()) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Load iClusion DB", options);
            System.exit(1);
        }

        LOGGER.info("Reading iClusion trial TSV from '{}'", iclusionDBTsv);
        List<IclusionTrial> trials = IclusionTrialFile.read(iclusionDBTsv);
        LOGGER.info(" Read {} trials", trials.size());

        IclusionDatabaseAccess dbWriter = IclusionDatabaseAccess.databaseAccess(cmd);

        dbWriter.writeIclusionDAO(trials);
        LOGGER.info("Written iclusion db to database");
    }

    private static Options createOptions() {
        Options options = new Options();
        options.addOption(ICLUSION_DB_TSV, true, "Path towards the iClusion db tsv.");

        IclusionDatabaseAccess.addDatabaseCmdLineArgs(options);
        return options;
    }
}