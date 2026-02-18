package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;

import java.nio.file.Path;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.finding.CurationApplier;
import com.hartwig.hmftools.datamodel.finding.CurationRecord;
import com.hartwig.hmftools.datamodel.finding.FindingRecord;
import com.hartwig.hmftools.datamodel.finding.FindingsJson;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CuratedFindingsApplication
{
    static class Config
    {
        public final String InputFindingJsonPath;
        public final String CurationJsonPath;
        public final String OutputFindingJsonPath;

        private static final String INPUT_FINDING_JSON_PATH_ARG = "input_findings";
        private static final String CURATION_JSON_PATH_ARG = "curations";
        private static final String OUTPUT_FINDING_JSON_PATH_ARG = "output_findings";

        public Config(final ConfigBuilder configBuilder)
        {
            this(configBuilder.getValue(INPUT_FINDING_JSON_PATH_ARG),
                    configBuilder.getValue(CURATION_JSON_PATH_ARG),
                    configBuilder.getValue(OUTPUT_FINDING_JSON_PATH_ARG));
        }

        public Config(final String inputFindingJsonPath, final String curationJsonPath, final String outputFindingJsonPath)
        {
            this.InputFindingJsonPath = inputFindingJsonPath;
            this.CurationJsonPath = curationJsonPath;
            this.OutputFindingJsonPath = outputFindingJsonPath;
        }

        public static void registerConfig(final ConfigBuilder configBuilder)
        {
            configBuilder.addPath(INPUT_FINDING_JSON_PATH_ARG, true, "Path to input findings JSON file");
            configBuilder.addPath(CURATION_JSON_PATH_ARG, true, "Path to curation JSON file");
            configBuilder.addConfigItem(OUTPUT_FINDING_JSON_PATH_ARG, true, "Path to output findings JSON file");

            addLoggingOptions(configBuilder);
        }
    }

    public static final String APP_NAME = "Finding";

    public static final Logger LOGGER = LogManager.getLogger(FindingApplication.class);

    public static void main(String[] args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        Config.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        Config config = new Config(configBuilder);
        new CuratedFindingsApplication(config).run();
    }

    private final Config config;

    private CuratedFindingsApplication(final Config config)
    {
        this.config = config;
    }

    private void run() throws Exception
    {
        LOGGER.info("Applying curations({}) to findings({}), output findings({})",
                config.CurationJsonPath, config.InputFindingJsonPath, config.OutputFindingJsonPath);

        final FindingRecord inputFindingRecord = new FindingsJson().read(Path.of(config.InputFindingJsonPath));
        final CurationRecord curationRecord = CurationRecord.jsonReadWriter().read(Path.of(config.CurationJsonPath));

        final FindingRecord outputFindingRecord = CurationApplier.applyCurations(inputFindingRecord, curationRecord);

        new FindingsJson().write(outputFindingRecord, Path.of(config.OutputFindingJsonPath));

        LOGGER.info("Done!");
    }
}
