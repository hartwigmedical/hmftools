package com.hartwig.hmftools.common.peach;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class PeachFactory {

    private static final Logger LOGGER = LogManager.getLogger(PeachFactory.class);

    private PeachFactory() {
    }

    @NotNull
    public static List<PeachGenotype> analyzePeach(@NotNull String peachGenotypeTsv) throws IOException {
        LOGGER.info("Loading pharmacogenetics data from {}", new File(peachGenotypeTsv).getParent());
        List<PeachGenotype> peachGenotypes = PeachGenotypeFile.read(peachGenotypeTsv);
        LOGGER.info(" Loaded {} reportable pharmacogenetics from {}", peachGenotypes.size(), peachGenotypeTsv);
        return peachGenotypes;
    }
}
