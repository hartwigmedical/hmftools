package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.BEDFileLoader;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadAmberData {

    private static final Logger LOGGER = LogManager.getLogger(LoadAmberData.class);

    private static final String SAMPLE = "sample";

    private static final String AMBER_DIR = "amber_dir";
    private static final String BED_FILE = "bed";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String amberPath = cmd.getOptionValue(AMBER_DIR);

        String amberFile = AmberBAFFile.generateAmberFilenameForReading(amberPath, tumorSample);
        LOGGER.info("Loading data from {}", amberFile);

        Multimap<Chromosome, AmberBAF> bafs = AmberBAFFile.read(amberFile);
        GenomePositionSelector<AmberBAF> selector = GenomePositionSelectorFactory.create(bafs);

        List<GenomeRegion> loci = Lists.newArrayList(BEDFileLoader.fromBedFile(cmd.getOptionValue(BED_FILE)).values());
        List<AmberBAF> lociAmberPoints = Lists.newArrayList();
        for (GenomeRegion locus : loci) {
            selector.select(GenomePositions.create(locus.chromosome(), locus.start())).ifPresent(lociAmberPoints::add);
        }

        persistToDatabase(cmd, tumorSample, lociAmberPoints);
        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample");
        options.addOption(AMBER_DIR, true, "Path to the amber directory");
        addDatabaseCmdLineArgs(options);
        options.addOption(BED_FILE, true, "Location of bed file");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    private static void persistToDatabase(@NotNull CommandLine cmd, @NotNull String tumorSample, @NotNull List<AmberBAF> amber)
            throws SQLException {
        try (DatabaseAccess dbAccess = databaseAccess(cmd)) {
            dbAccess.writeAmberBAF(tumorSample, amber);
        }
    }
}
