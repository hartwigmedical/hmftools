package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.hla.HlaFiles;
import com.hartwig.hmftools.common.hla.HlaType;
import com.hartwig.hmftools.common.hla.HlaTypeDetails;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadLilacData {

    private static final Logger LOGGER = LogManager.getLogger(LoadLilacData.class);

    private static final String SAMPLE = "sample";
    private static final String LILAC_DIR = "lilac_dir";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        String sample = cmd.getOptionValue(SAMPLE);
        String lilacDir = cmd.getOptionValue(LILAC_DIR);

        LOGGER.info("Reading data from {}", lilacDir);
        final String lilacFile = lilacDir + "/" + sample + ".lilac.txt";
        final String qcFile = lilacDir + "/" + sample + ".lilac.qc.txt";
        final HlaType type = HlaFiles.type(lilacFile, qcFile);
        final List<HlaTypeDetails> details = HlaFiles.typeDetails(lilacFile);

        LOGGER.info("Persisting lilac data to db for {}", sample);
        dbAccess.writeHla(sample, type, details);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(LILAC_DIR, true, "Path to the lilac directory.");
        addDatabaseCmdLineArgs(options);
        return options;
    }

}
