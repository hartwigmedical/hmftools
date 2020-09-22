package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.pharmacogenetics.PGXCalls;
import com.hartwig.hmftools.common.pharmacogenetics.PGXCallsFile;
import com.hartwig.hmftools.common.pharmacogenetics.PGXGenotype;
import com.hartwig.hmftools.common.pharmacogenetics.PGXGenotypeFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadPgxData {
    private static final Logger LOGGER = LogManager.getLogger(LoadPgxData.class);

    private static final String SAMPLE = "sample";
    private static final String PGX_CALLS_TXT = "pgx_calls_txt";
    private static final String PGX_GENOTYPE_TXT = "pgx_genotype_txt";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);

        String pgxCallsFileName = cmd.getOptionValue(PGX_CALLS_TXT);
        String pgxGenotypeFileName = cmd.getOptionValue(PGX_GENOTYPE_TXT);

        if (Utils.anyNull(sample, pgxCallsFileName, pgxGenotypeFileName)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load metrics data", options);
        } else {
            DatabaseAccess dbWriter = databaseAccess(cmd);

            LOGGER.info("Reading pgx calls file {}", pgxCallsFileName);
            List<PGXCalls> pgxCalls = PGXCallsFile.read(pgxCallsFileName);
            LOGGER.info(" Read {} pgx calls", pgxCalls.size());

            LOGGER.info("Reading pgx genotype file {}", pgxGenotypeFileName);
            List<PGXGenotype> pgxGenotype = PGXGenotypeFile.read(pgxGenotypeFileName);
            LOGGER.info(" Read {} pgx genotypes", pgxGenotype.size());

            LOGGER.info("Writing pgx into database for {}", sample);
            dbWriter.writePGX(sample, pgxGenotype, pgxCalls);
            LOGGER.info("Complete");
        }
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample for which we are going to load the metrics");

        options.addOption(PGX_CALLS_TXT, true, "Path towards the pgx calls txt file");
        options.addOption(PGX_GENOTYPE_TXT, true, "Path towards the pgx genotype txt file");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
