package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadGermlineDriverCatalog {

    private static final Logger LOGGER = LogManager.getLogger(LoadGermlineDriverCatalog.class);

    private static final String SAMPLE = "sample";
    private static final String PURPLE_DIR = "purple_dir";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String purpleDir = cmd.getOptionValue(PURPLE_DIR);

        String germlineDriverCatalogFileName = DriverCatalogFile.generateGermlineFilename(purpleDir, tumorSample);
        List<DriverCatalog> germlineDriverCatalog = DriverCatalogFile.read(germlineDriverCatalogFileName);
        LOGGER.info("Read {} drivers from {}", germlineDriverCatalog.size(), germlineDriverCatalogFileName);

        dbAccess.writeGermlineDriverCatalog(tumorSample, germlineDriverCatalog);
        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(PURPLE_DIR, true, "Path to the purple directory.");
        addDatabaseCmdLineArgs(options);
        return options;
    }
}
