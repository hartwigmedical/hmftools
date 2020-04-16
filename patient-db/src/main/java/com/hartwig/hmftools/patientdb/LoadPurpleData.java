package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadPurpleData {

    private static final Logger LOGGER = LogManager.getLogger(LoadPurpleData.class);

    private static final String SAMPLE = "sample";

    private static final String PURPLE_DIR = "purple_dir";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String purplePath = cmd.getOptionValue(PURPLE_DIR);

        LOGGER.info("Reading data from {}", purplePath);
        PurityContext purityContext = FittedPurityFile.read(purplePath, tumorSample);
        List<GeneCopyNumber> geneCopyNumbers =
                GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilenameForReading(purplePath, tumorSample));
        List<PurpleCopyNumber> copyNumbers =
                PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(purplePath, tumorSample));
        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(DriverCatalogFile.generateFilename(purplePath, tumorSample));

        PurpleQC purpleQC = PurpleQCFile.read(PurpleQCFile.generateFilename(purplePath, tumorSample));
        List<FittedPurity> bestFitPerPurity = FittedPurityRangeFile.readBestFitPerPurity(purplePath, tumorSample);

        String germlineCopyNumberFilename = PurpleCopyNumberFile.generateGermlineFilenameForReading(purplePath, tumorSample);
        List<PurpleCopyNumber> germlineCopyNumbers = new File(germlineCopyNumberFilename).exists()
                ? PurpleCopyNumberFile.read(germlineCopyNumberFilename)
                : Lists.newArrayList();

        LOGGER.info("Persisting to db");
        persistToDatabase(dbAccess,
                tumorSample,
                bestFitPerPurity,
                copyNumbers,
                germlineCopyNumbers,
                purityContext,
                purpleQC,
                geneCopyNumbers,
                driverCatalog);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(PURPLE_DIR, true, "Path to the purple directory.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException {
        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }

    public static void persistToDatabase(@NotNull DatabaseAccess dbAccess, @NotNull String tumorSample,
            @NotNull List<FittedPurity> bestFitPerPurity, @NotNull List<PurpleCopyNumber> copyNumbers,
            @NotNull List<PurpleCopyNumber> germlineDeletions, @NotNull PurityContext purityContext, @NotNull PurpleQC qcChecks,
            @NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull List<DriverCatalog> driverCatalog) {
        dbAccess.writePurity(tumorSample, purityContext, qcChecks);
        dbAccess.writeBestFitPerPurity(tumorSample, bestFitPerPurity);
        dbAccess.writeCopynumbers(tumorSample, copyNumbers);
        dbAccess.writeGermlineCopynumbers(tumorSample, germlineDeletions);
        dbAccess.writeGeneCopynumberRegions(tumorSample, geneCopyNumbers);
        dbAccess.writeDriverCatalog(tumorSample, driverCatalog);
    }
}
