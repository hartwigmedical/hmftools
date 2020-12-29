package com.hartwig.hmftools.serve.extraction.codon.tools;


import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CodonChecker {
    private static final Logger LOGGER = LogManager.getLogger(CodonChecker.class);

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE checker on version {}");

        LOGGER.info("Checking codons");

        LOGGER.info("Done!");

    }

}
