package com.hartwig.hmftools.orange;

import java.io.IOException;

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

    public static void main(String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        OrangeConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        OrangeConfig config = OrangeConfig.createConfig(configBuilder);
        new OrangeApplication(config).run();
    }

    private final OrangeConfig config;

    private OrangeApplication(final OrangeConfig config)
    {
        this.config = config;
    }

    private void run() throws IOException
    {
        LOGGER.info("Generating ORANGE report data");

        OrangeAlgo algo = OrangeAlgo.fromConfig(config);
        OrangeRecord report = algo.run(config);

        ReportWriter writer = ReportWriterFactory.createToDiskWriter(config);
        writer.write(report);

        LOGGER.info("Done!");
    }
}