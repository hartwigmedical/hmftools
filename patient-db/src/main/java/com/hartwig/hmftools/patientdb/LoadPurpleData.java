package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
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

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final DatabaseAccess dbAccess = databaseAccess(cmd);

        final String tumorSample = cmd.getOptionValue(SAMPLE);
        final String purplePath = cmd.getOptionValue(PURPLE_DIR);

        LOGGER.info("Reading data from {}", purplePath);
        final PurpleQC purpleQC = PurpleQCFile.read(PurpleQCFile.generateFilename(purplePath, tumorSample));
        final PurityContext purityContext = FittedPurityFile.read(purplePath, tumorSample);
        final List<FittedPurity> bestFitPerPurity = FittedPurityRangeFile.read(purplePath, tumorSample);
        final List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilename(purplePath, tumorSample));
        final List<PurpleCopyNumber> copyNumbers =
                PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilename(purplePath, tumorSample));
        final List<PurpleCopyNumber> germlineDeletions =
                PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateGermlineFilename(purplePath, tumorSample));
        final List<FittedRegion> enrichedFittedRegions = FittedRegionFile.read(FittedRegionFile.generateFilename(purplePath, tumorSample));

        LOGGER.info("Persisting to db");
        persistToDatabase(dbAccess,
                tumorSample,
                bestFitPerPurity,
                copyNumbers,
                germlineDeletions,
                enrichedFittedRegions,
                purityContext,
                purpleQC,
                geneCopyNumbers);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(PURPLE_DIR, true, "Path to the purple directory.");
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

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }

    public static void persistToDatabase(final DatabaseAccess dbAccess, final String tumorSample, final List<FittedPurity> bestFitPerPurity,
            final List<PurpleCopyNumber> copyNumbers, final List<PurpleCopyNumber> germlineDeletions,
            final List<FittedRegion> enrichedFittedRegions, final PurityContext purityContext, final PurpleQC qcChecks,
            final List<GeneCopyNumber> geneCopyNumbers) {
        dbAccess.writePurity(tumorSample, purityContext, qcChecks);
        dbAccess.writeBestFitPerPurity(tumorSample, bestFitPerPurity);
        dbAccess.writeCopynumbers(tumorSample, copyNumbers);
        dbAccess.writeGermlineCopynumbers(tumorSample, germlineDeletions);
        dbAccess.writeCopynumberRegions(tumorSample, enrichedFittedRegions);
        dbAccess.writeGeneCopynumberRegions(tumorSample, geneCopyNumbers);
    }

}
