package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import java.nio.file.Path;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingsJson;
import com.hartwig.hmftools.finding.datamodel.JsonReadWriter;
import com.hartwig.hmftools.finding.datamodel.VisualisationFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VisualisationFileFinderApp
{
    public static class Config
    {
        public final Path FindingJsonPath;
        public final Path OutputJsonPath;

        private static final String FINDING_JSON_PATH_ARG = "finding_json";
        private static final String OUTPUT_JSON_PATH_ARG = "output_json";

        public Config(final ConfigBuilder configBuilder)
        {
            this(configBuilder.getValue(FINDING_JSON_PATH_ARG),  configBuilder.getValue(OUTPUT_JSON_PATH_ARG));
        }

        public Config(final String findingJsonPath, final String outputJsonPath)
        {
            this.FindingJsonPath = Path.of(findingJsonPath);
            this.OutputJsonPath = Path.of(outputJsonPath);
        }

        public static void registerConfig(final ConfigBuilder configBuilder)
        {
            configBuilder.addConfigItem(FINDING_JSON_PATH_ARG, true, "Path to input finding JSON file");
            configBuilder.addConfigItem(OUTPUT_JSON_PATH_ARG, true, "Path to output visualisation file list JSON file");
            addLoggingOptions(configBuilder);
        }
    }

    public static final String APP_NAME = "Finding";

    public static final Logger LOGGER = LogManager.getLogger(VisualisationFileFinderApp.class);

    public static void main(String[] args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        Config.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        Config config = new Config(configBuilder);
        new VisualisationFileFinderApp(config).run();
    }

    private final Config config;

    private VisualisationFileFinderApp(final Config config)
    {
        this.config = config;
    }

    private void run() throws Exception
    {
        LOGGER.info("Visualisation file finder");

        FindingRecord findingRecord = new FindingsJson().read(config.FindingJsonPath);
        List<String> visualisationFiles = VisualisationFileFinder.find(findingRecord).stream()
                .map(VisualisationFile::fileName)
                .toList();

        JsonReadWriter.listOf(String.class).write(visualisationFiles, config.OutputJsonPath);

        LOGGER.info("Done!");
    }
}
