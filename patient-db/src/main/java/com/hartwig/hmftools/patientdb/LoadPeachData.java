package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.peach.PeachCalls;
import com.hartwig.hmftools.common.peach.PeachCallsFile;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadPeachData {

    private static final Logger LOGGER = LogManager.getLogger(LoadPeachData.class);

    private static final String SAMPLE = "sample";
    private static final String PEACH_CALLS_TXT = "peach_calls_txt";
    private static final String PEACH_GENOTYPE_TXT = "peach_genotype_txt";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);

        String peachCallsFileName = cmd.getOptionValue(PEACH_CALLS_TXT);
        String peachGenotypeFileName = cmd.getOptionValue(PEACH_GENOTYPE_TXT);

        if (Utils.anyNull(sample, peachCallsFileName, peachGenotypeFileName)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load PGX Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Reading PEACH calls file {}", peachCallsFileName);
        List<PeachCalls> peachCalls = PeachCallsFile.read(peachCallsFileName);
        LOGGER.info(" Read {} PEACH calls", peachCalls.size());

        LOGGER.info("Reading PEACH genotype file {}", peachGenotypeFileName);
        List<PeachGenotype> peachGenotype = PeachGenotypeFile.read(peachGenotypeFileName);
        LOGGER.info(" Read {} PEACH genotypes", peachGenotype.size());

        LOGGER.info("Writing PEACH into database for {}", sample);
        dbWriter.writePeach(sample, peachGenotype, peachCalls);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample for which we are going to load the metrics");

        options.addOption(PEACH_CALLS_TXT, true, "Path towards the PEACH calls txt file");
        options.addOption(PEACH_GENOTYPE_TXT, true, "Path towards the PEACH genotype txt file");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
