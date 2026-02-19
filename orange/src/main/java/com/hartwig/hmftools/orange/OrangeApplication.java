package com.hartwig.hmftools.orange;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.report.ReportWriter;
import com.hartwig.hmftools.orange.report.ReportWriterFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class OrangeApplication
{
    public static final String APP_NAME = "Orange";

    public static final Logger LOGGER = LogManager.getLogger(OrangeApplication.class);

    public static final String VERSION = OrangeApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        OrangeConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        OrangeConfig config = new OrangeConfig(configBuilder);
        new OrangeApplication(config).run();
    }

    private final OrangeConfig mConfig;

    private OrangeApplication(final OrangeConfig config)
    {
        mConfig = config;
    }

    private void run() throws Exception
    {
        LOGGER.info("Generating ORANGE report data");

        OrangeAlgo algo = OrangeAlgo.fromConfig(mConfig);
        OrangeRecord orangeRecord = algo.run(mConfig);

        ReportWriter writer = ReportWriterFactory.createToDiskWriter(mConfig);
        writer.write(orangeRecord);

        LOGGER.info("Done!");
    }
}