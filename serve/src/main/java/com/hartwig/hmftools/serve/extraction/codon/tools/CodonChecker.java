package com.hartwig.hmftools.serve.extraction.codon.tools;

import java.io.IOException;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class CodonChecker {
    private static final Logger LOGGER = LogManager.getLogger(CodonChecker.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE codon checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String knownCodonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownCodons.SERVE.37.tsv";

        LOGGER.info("Checking codons");

        LOGGER.info("Done!");

    }

}
