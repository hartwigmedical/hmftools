package com.hartwig.hmftools.amber.blacklist;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class AmberTrainingConfig
{
    public final List<String> SampleIds;
    public final String AmberDir;
    public final String OutputFile;

    private static final String OUTPUT_FILE = "output_file";

    public AmberTrainingConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = Lists.newArrayList();
        loadSampleIds(configBuilder);

        AmberDir = configBuilder.getValue(AMBER_DIR_CFG, "");
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
    }

    private void loadSampleIds(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            return;
        }

        try
        {
            SampleIds.addAll(Files.readAllLines(Paths.get(configBuilder.getValue(SAMPLE_ID_FILE))));
        }
        catch(IOException e)
        {
            AMB_LOGGER.error("failed to load sample IDs: {}", e.toString());
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_ID_FILE, true, "CSV with SampleId, optional: WgsSampleId,Gender");
        configBuilder.addPath(AMBER_DIR_CFG, false, AMBER_DIR_DESC);
        configBuilder.addRequiredConfigItem(OUTPUT_FILE, "Output suspicious points file");
        addLoggingOptions(configBuilder);
    }
}
