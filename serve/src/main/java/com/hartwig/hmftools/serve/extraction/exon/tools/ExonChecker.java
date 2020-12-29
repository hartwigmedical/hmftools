package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ExonChecker {
    private static final Logger LOGGER = LogManager.getLogger(ExonChecker.class);

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE checker on version {}");

        LOGGER.info("Checking exons");

        LOGGER.info("Done!");

    }

}
