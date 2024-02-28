package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaPrediction;
import com.hartwig.hmftools.common.cuppa.interpretation.CuppaPredictionFactory;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadCuppa
{
    private static final String CUPPA_RESULTS_CSV = "cuppa_results_csv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(CUPPA_RESULTS_CSV, true, "Path to the *.cup.data.csv file");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);
        String cuppaResultsCsv = configBuilder.getValue(CUPPA_RESULTS_CSV);

        try(DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
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
        catch(Exception e)
        {
            LOGGER.error("Failed to load CUPPA", e);
            System.exit(1);
        }
    }
}