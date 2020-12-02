package com.hartwig.hmftools.serve.sources.vicc.check;

import java.util.Set;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GeneChecker {

    private static final Logger LOGGER = LogManager.getLogger(GeneChecker.class);

    @NotNull
    private final Set<String> allValidGenes;

    @NotNull
    public static GeneChecker buildForHG19() {
        return new GeneChecker(HmfGenePanelSupplier.allGenesMap37().keySet());
    }

    public GeneChecker(@NotNull final Set<String> allValidGenes) {
        this.allValidGenes = allValidGenes;
    }

    public boolean isValidGene(@Nullable String gene) {
        if (allValidGenes.contains(gene)) {
            return true;
        } else {
            if (gene != null) {
                LOGGER.warn("Gene '{}' is not defined in the exome used!", gene);
            }
            return false;
        }
    }
}
