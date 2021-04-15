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
        if (allValidGenes.contains(gene)) {
            return true;
        } else {
            if (gene != null) {
                LOGGER.warn("Gene '{}' is not present in the full gene list used!", gene);
            }
            return false;
        }
    }
}
