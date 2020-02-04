package com.hartwig.hmftools.knowledgebasegenerator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class KnowledgebaseGeneratorApplication {

    private static final Logger LOGGER = LogManager.getLogger(KnowledgebaseGeneratorApplication.class);

    private static final String ICLUSION_TRIAL_TSV = "iclusion_trial_tsv";
    private static final String VICC_JSON = "vicc_json";

    private static final String VERSION = KnowledgebaseGeneratorApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) {
        LOGGER.info("Running Knowledgebase Generator v{}", VERSION);

    }
}
