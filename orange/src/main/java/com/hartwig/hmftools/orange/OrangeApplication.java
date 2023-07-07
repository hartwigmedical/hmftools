package com.hartwig.hmftools.orange;

import java.io.IOException;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.report.ReportWriter;
import com.hartwig.hmftools.orange.report.ReportWriterFactory;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class OrangeApplication {

    private static final Logger LOGGER = LogManager.getLogger(OrangeApplication.class);

    private static final String APPLICATION = "ORANGE";
    public static final String VERSION = OrangeApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException, ParseException {
        LOGGER.info("Running {} v{}", APPLICATION, VERSION);

        ConfigBuilder configBuilder = new ConfigBuilder();
        OrangeConfig.registerConfig(configBuilder);

        if (!configBuilder.parseCommandLine(args)) {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        OrangeConfig config = OrangeConfig.createConfig(configBuilder);
        new OrangeApplication(config).run();
    }

    @NotNull
    private final OrangeConfig config;

    private OrangeApplication(@NotNull final OrangeConfig config) {
        this.config = config;
    }

    private void run() throws IOException {
        LOGGER.info("Generating ORANGE report data");
        OrangeAlgo algo = OrangeAlgo.fromConfig(config);
        OrangeRecord report = algo.run(config);

        ReportWriter writer = ReportWriterFactory.createToDiskWriter(config);
        writer.write(report);

        LOGGER.info("Done!");
    }
}