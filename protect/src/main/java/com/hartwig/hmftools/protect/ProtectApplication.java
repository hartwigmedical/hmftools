package com.hartwig.hmftools.protect;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProtectApplication {

    private static final Logger LOGGER = LogManager.getLogger(ProtectApplication.class);

    public static void main(@NotNull String[] args) {
        Options options = ProtectConfig.createOptions();

        ProtectConfig config = null;
        try {
            config = ProtectConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PROTECT", options);
            System.exit(1);
        }

        LOGGER.info("Running PROTECT with {}", config);
    }
}
