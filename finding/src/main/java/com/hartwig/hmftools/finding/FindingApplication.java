package com.hartwig.hmftools.finding;

import java.nio.file.Path;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.finding.FindingRecord;
import com.hartwig.hmftools.datamodel.finding.FindingsJson;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FindingApplication
{
    public static final String APP_NAME = "Finding";

    public static final Logger LOGGER = LogManager.getLogger(FindingApplication.class);

    public static void main(String[] args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        FindingApplicationConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        FindingApplicationConfig config = new FindingApplicationConfig(configBuilder);
        new FindingApplication(config).run();
    }

    private final FindingApplicationConfig config;

    private FindingApplication(final FindingApplicationConfig config)
    {
        this.config = config;
    }

    private void run() throws Exception
    {
        LOGGER.info("Generating finding json");

        OrangeRecord orangeRecord = OrangeJson.getInstance().read(config.OrangeJsonPath);

        FindingRecord findingRecord = FindingRecordFactory.fromOrangeRecord(orangeRecord,
                config.ClinicalTranscriptsPath != null ? Path.of(config.ClinicalTranscriptsPath) : null,
                Path.of(config.DriverGenePath));

        FindingsJson.getInstance().write(findingRecord, Path.of(config.FindingJsonPath));
        LOGGER.info("Done!");
    }
}
