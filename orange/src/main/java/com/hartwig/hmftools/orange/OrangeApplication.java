package com.hartwig.hmftools.orange;

import java.io.IOException;

import com.hartwig.hmftools.orange.algo.OrangeAlgo;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportWriter;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class OrangeApplication {

    private static final Logger LOGGER = LogManager.getLogger(OrangeApplication.class);

    public static final String VERSION = OrangeApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running ORANGE v{}", VERSION);

        Options options = OrangeConfig.createOptions();

        OrangeConfig config = null;
        try {
            config = OrangeConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("ORANGE", options);
            System.exit(1);
        }

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
        OrangeReport report = algo.run(config);

        LOGGER.info("Writing report");
        ReportWriter writer = new ReportWriter(true, config.outputDir(), config.reportGermline());
        writer.write(report);

        LOGGER.info("Done!");
    }
}