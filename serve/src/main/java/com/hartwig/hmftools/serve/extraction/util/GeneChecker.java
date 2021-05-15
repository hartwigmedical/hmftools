package com.hartwig.hmftools.serve.extraction.util;

import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneChecker {

    private static final Logger LOGGER = LogManager.getLogger(GeneChecker.class);

    @NotNull
    private final Set<String> allValidGenes;

    public GeneChecker(@NotNull final Set<String> allValidGenes) {
        this.allValidGenes = allValidGenes;
    }

    public boolean isValidGene(@Nullable String gene) {
        if (geneExistsInAllValidGenes(gene)) {
            return true;
        } else {
            if (gene != null) {
                LOGGER.warn("Gene '{}' is not considered a valid gene!", gene);
            }
            return false;
        }
    }

    public boolean geneExistsInAllValidGenes(@Nullable String gene) {
        return allValidGenes.contains(gene);
    }
}
