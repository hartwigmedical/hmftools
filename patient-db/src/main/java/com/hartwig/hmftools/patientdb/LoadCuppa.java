package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadCuppa
{
    private static final String CUPPA_VIS_DATA_TSV = "cuppa_vis_data_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        String sample = cmd.getOptionValue(SAMPLE);

        logVersion();

        String cuppaVisDataTsv = cmd.getOptionValue(CUPPA_VIS_DATA_TSV);

        if(CommonUtils.anyNull(sample, cuppaVisDataTsv))
        {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load CUPPA Data", options);
            System.exit(1);
        }

        LOGGER.info("Loading CUPPA from {}", new File(cuppaVisDataTsv).getParent());
        CuppaPredictions cuppaPredictions = CuppaPredictions.fromTsv(cuppaVisDataTsv);
        LOGGER.info(" Loaded {} entries from {}", cuppaPredictions.size(), cuppaVisDataTsv);

//        List<CuppaPrediction> predictions = CuppaPredictionFactory.create(cuppaEntries);
//        CuppaPrediction best = predictions.get(0);
//        LOGGER.info(" Predicted cancer type '{}' with likelihood {}", best.cancerType(), best.likelihood());

        LOGGER.info("Writing CUPPA into database for {}", sample);
        DatabaseAccess dbWriter = databaseAccess(cmd);

        for(CuppaPredictionEntry cuppaPredictionEntry: cuppaPredictions.PredictionEntries)
        {
            dbWriter.writeCuppa(
                    cuppaPredictionEntry.SampleId,
                    cuppaPredictionEntry.DataType.toString(),
                    cuppaPredictionEntry.ClfName.toString(),
                    cuppaPredictionEntry.CancerType,
                    cuppaPredictionEntry.DataValue,
                    cuppaPredictionEntry.Rank,
                    cuppaPredictionEntry.RankGroup
            );
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample name");
        options.addOption(CUPPA_VIS_DATA_TSV, true, "Path to the CUPPA vis data file");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}