package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
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

        String sample = cmd.getOptionValue(SAMPLE);
        String purpleDir = cmd.getOptionValue(PURPLE_DIR);

        LOGGER.info("Reading data from {}", purpleDir);
        PurityContext purityContext = PurityContextFile.read(purpleDir, sample);
        List<FittedPurity> bestFitPerPurity = FittedPurityRangeFile.readBestFitPerPurity(purpleDir, sample);

        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilenameForReading(purpleDir, sample));
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(purpleDir, sample));
        List<GermlineDeletion> germlineDeletions = GermlineDeletion.read(GermlineDeletion.generateFilename(purpleDir, sample));

        List<DriverCatalog> somaticDriverCatalog = DriverCatalogFile.read(DriverCatalogFile.generateSomaticFilename(purpleDir, sample));

        final String germlineDriverFile = DriverCatalogFile.generateGermlineFilename(purpleDir, sample);

        List<DriverCatalog> germlineDriverCatalog =
                Files.exists(Paths.get(germlineDriverFile)) ? DriverCatalogFile.read(germlineDriverFile) : Lists.newArrayList();

        LOGGER.info("Persisting purple data to db for {}", sample);

        dbAccess.writePurity(sample, purityContext, purityContext.qc());
        dbAccess.writeBestFitPerPurity(sample, bestFitPerPurity);
        dbAccess.writeCopynumbers(sample, copyNumbers);
        dbAccess.writeGeneCopyNumbers(sample, geneCopyNumbers);
        dbAccess.writeGermlineDeletions(sample, germlineDeletions);
        dbAccess.writePurpleDriverCatalog(sample, somaticDriverCatalog, germlineDriverCatalog);

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
