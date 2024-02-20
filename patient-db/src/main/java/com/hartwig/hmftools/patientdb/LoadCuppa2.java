package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.*;
import org.jetbrains.annotations.NotNull;

public class LoadCuppa2
{
    private static final String CUPPA_VIS_DATA_TSV = "cuppa_vis_data_tsv";

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();

        options.addOption(CUPPA_VIS_DATA_TSV, true, "Path to the CUPPA vis data file");
        options.addOption(SAMPLE, true, "Sample for which we are going to load the CUPPA results");

        addDatabaseCmdLineArgs(options);

        return options;
    }

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        String cuppaVisDataTsv = cmd.getOptionValue(CUPPA_VIS_DATA_TSV);

        if(CommonUtils.anyNull(cuppaVisDataTsv))
        {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load CUPPA Data", options);
            System.exit(1);
        }

        logVersion();

        try
        {
            LOGGER.info("Loading CUPPA from {}", new File(cuppaVisDataTsv).getParent());
            CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(cuppaVisDataTsv);
            String sample = cmd.getOptionValue(SAMPLE);
            LOGGER.info("Loaded {} entries from {} for sample {}", cuppaPredictions.size(), cuppaVisDataTsv, sample);

            int TOP_N_PROBS = 3;
            LOGGER.info("Writing top {} probabilities from all classifiers to database", TOP_N_PROBS);

            try (DatabaseAccess dbWriter = databaseAccess(cmd))
            {
                dbWriter.writeCuppa2(sample, cuppaPredictions, TOP_N_PROBS);
                LOGGER.info("Complete");
            }
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load CUPPA data", e);
            System.exit(1);
        }
    }
}