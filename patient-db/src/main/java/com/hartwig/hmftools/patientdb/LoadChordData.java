package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.chord.ChordDataFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadChordData
{
    private static final String PREDICTION_FILE = "prediction_file";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(PREDICTION_FILE, true, "Path towards the chord prediction file.");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);
        String predictionFile = configBuilder.getValue(PREDICTION_FILE);

        try (DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("Extracting and writing chord for {}", predictionFile);
            ChordData chordData = ChordDataFile.read(predictionFile);
            dbWriter.writeChord(sample, chordData);

            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load Chord data", e);
            System.exit(1);
        }
    }
}