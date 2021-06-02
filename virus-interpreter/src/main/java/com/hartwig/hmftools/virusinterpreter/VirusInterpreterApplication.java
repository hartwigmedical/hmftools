package com.hartwig.hmftools.virusinterpreter;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VirusInterpreterApplication {

    private static final Logger LOGGER = LogManager.getLogger(VirusInterpreterApplication.class);

    private static final String VERSION = VirusInterpreterApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) {
        LOGGER.info("Running Virus Interpreter v{}", VERSION);

        Options options = VirusInterpreterConfig.createOptions();

        VirusInterpreterConfig config = null;
        try {
            config = VirusInterpreterConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("VirusInterpreter", options);
            System.exit(1);
        }

        LOGGER.info(config);
    }
}
