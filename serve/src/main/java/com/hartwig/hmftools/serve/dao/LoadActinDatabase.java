package com.hartwig.hmftools.serve.dao;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinFileReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadActinDatabase {

    private static final Logger LOGGER = LogManager.getLogger(LoadActinDatabase.class);

    private static final String ACTIN_DB_TSV = "actin_db_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String actinDBTsv = cmd.getOptionValue(ACTIN_DB_TSV);

        if (Utils.anyNull(actinDBTsv) || !new File(actinDBTsv).exists()) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Load ACTIN DB", options);
            System.exit(1);
        }

        LOGGER.info("Reading ACTIN trial TSV from '{}'", actinDBTsv);
        List<ActinEntry> trials = ActinFileReader.read(actinDBTsv);
        LOGGER.info(" Read {} trials", trials.size());

        KnowledgebaseDatabaseAccess dbWriter = KnowledgebaseDatabaseAccess.databaseAccess(cmd);

        dbWriter.writeActinDAO(trials);
        LOGGER.info("Written ACTN db to database");
    }

    private static Options createOptions() {
        Options options = new Options();
        options.addOption(ACTIN_DB_TSV, true, "Path towards the ACTIN db tsv.");

        KnowledgebaseDatabaseAccess.addDatabaseCmdLineArgs(options);
        return options;
    }
}