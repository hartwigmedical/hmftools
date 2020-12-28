package com.hartwig.hmftools.serve;

import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ServeCheckerApplication {

    private static final Logger LOGGER = LogManager.getLogger(ServeApplication.class);
    private static final String VERSION = ServeApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE checker on version {}", VERSION);

        LOGGER.info("Checking codons");

        LOGGER.info("Checking exons");

        LOGGER.info("Done!");

    }
}
