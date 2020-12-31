package com.hartwig.hmftools.serve.extraction.exon.tools;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExonFile;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

public class ExonChecker {
    private static final Logger LOGGER = LogManager.getLogger(ExonChecker.class);
    private static final boolean LOG_DEBUG = true;

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running SERVE exon checker");

        if (LOG_DEBUG) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String knownExonsTsv = System.getProperty("user.home") + "/hmf/tmp/serve/KnownExons.SERVE.37.tsv";
        List<KnownExon> exons = KnownExonFile.read(knownExonsTsv);

        LOGGER.info("The size of the file is {}", exons.size());

        for (KnownExon exon: exons) {
            LOGGER.info(exon);
        }
        LOGGER.info("Checking exons");

        LOGGER.info("Done!");

    }

}
