package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
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

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String purplePath = cmd.getOptionValue(PURPLE_DIR);

        LOGGER.info("Reading data from {}", purplePath);
        PurityContext purityContext = PurityContextFile.read(purplePath, tumorSample);
        List<GeneCopyNumber> geneCopyNumbers =
                GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilenameForReading(purplePath, tumorSample));
        List<PurpleCopyNumber> copyNumbers =
                PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(purplePath, tumorSample));
        List<DriverCatalog> somaticDriverCatalog =
                DriverCatalogFile.read(DriverCatalogFile.generateSomaticFilenameForReading(purplePath, tumorSample));

        String germlineDriverCatalogFileName = DriverCatalogFile.generateGermlineFilename(purplePath, tumorSample);
        List<DriverCatalog> germlineDriverCatalog = new File(germlineDriverCatalogFileName).exists()
                ? DriverCatalogFile.read(germlineDriverCatalogFileName)
                : Collections.emptyList();

        List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(somaticDriverCatalog);
        driverCatalog.addAll(germlineDriverCatalog);

        PurpleQC purpleQC = purityContext.qc();
        List<FittedPurity> bestFitPerPurity = FittedPurityRangeFile.readBestFitPerPurity(purplePath, tumorSample);

        String germlineCopyNumberFilename = PurpleCopyNumberFile.generateGermlineFilenameForReading(purplePath, tumorSample);
        List<PurpleCopyNumber> germlineCopyNumbers = new File(germlineCopyNumberFilename).exists()
                ? PurpleCopyNumberFile.read(germlineCopyNumberFilename)
                : Lists.newArrayList();

        LOGGER.info("Persisting purple data to db for {}", tumorSample);
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
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(PURPLE_DIR, true, "Path to the purple directory.");
        addDatabaseCmdLineArgs(options);
        return options;
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
