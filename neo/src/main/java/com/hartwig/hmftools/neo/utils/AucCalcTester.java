package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.stats.AucCalc;
import com.hartwig.hmftools.common.stats.AucData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.NotNull;

public class AucCalcTester
{
    private static final String DATA_FILE = "data_file";
    private static final String IS_PERCENTILE = "percentile";

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        configBuilder.addPath(DATA_FILE, true, "Input filename");
        configBuilder.addFlag(IS_PERCENTILE, "In percentile or rank form");
        configBuilder.addConfigItem(OUTPUT_DIR, true, "Output directory");
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        String dataFile = configBuilder.getValue(DATA_FILE);
        boolean isPercentile = configBuilder.hasFlag(IS_PERCENTILE);

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(dataFile));

            lines.remove(0);

            List<AucData> aucDataList = Lists.newArrayList();

            for(String line : lines)
            {
                String[] items = line.split(",");
                boolean isPositive = Boolean.parseBoolean(items[0]);
                double value = Double.parseDouble(items[1]);

                aucDataList.add(new AucData(isPositive, value, isPercentile));
            }

            NE_LOGGER.info("loaded {} AUC items from file({})", aucDataList.size(), dataFile);

            double auc = isPercentile ?
                    AucCalc.calcPercentilesAuc(aucDataList, Level.DEBUG) : AucCalc.calcScoresAuc(aucDataList, Level.DEBUG);

            NE_LOGGER.info(String.format("AUC from %s = %.4f",
                    isPercentile ? "percentiles" : "values", auc));
        }
        catch(IOException e)
        {
            NE_LOGGER.error("file load failed: {}", e.toString());
        }
    }
}
