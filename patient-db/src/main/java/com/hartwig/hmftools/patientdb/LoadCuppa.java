package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaPrediction;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaPredictionFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadCuppa {

    private static final Logger LOGGER = LogManager.getLogger(LoadCuppa.class);

    private static final String SAMPLE = "sample";
    private static final String CUPPA_RESULTS_CSV = "cuppa_results_csv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String cuppaResultsCsv = cmd.getOptionValue(CUPPA_RESULTS_CSV);

        if (Utils.anyNull(sample, cuppaResultsCsv)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load CUPPA Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Loading CUPPA from {}", new File(cuppaResultsCsv).getParent());
        List<CuppaDataFile> cuppaEntries = CuppaDataFile.read(cuppaResultsCsv);
        LOGGER.info(" Loaded {} entries from {}", cuppaEntries.size(), cuppaResultsCsv);

        List<CuppaPrediction> predictions = CuppaPredictionFactory.create(cuppaEntries);
        CuppaPrediction best = predictions.get(0);
        LOGGER.info(" Predicted cancer type '{}' with likelihood {}", best.cancerType(), best.likelihood());

        LOGGER.info("Writing CUPPA into database for {}", sample);
        dbWriter.writeCuppa(sample, best.cancerType(), best.likelihood());

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample for which we are going to load the CUPPA results");
        options.addOption(CUPPA_RESULTS_CSV, true, "Path towards the CUPPA conclusion txt file");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}