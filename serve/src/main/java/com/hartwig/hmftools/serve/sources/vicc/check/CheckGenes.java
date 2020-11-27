package com.hartwig.hmftools.serve.sources.vicc.check;

import java.util.Set;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CheckGenes {

    private static final Logger LOGGER = LogManager.getLogger(CheckGenes.class);

    private static final Set<String> GENES = Sets.newHashSet("IGH", "IGK", "IGL");

    public static void checkGensInPanel(@NotNull String gene, @NotNull String event) {
        if (!GENES.contains(gene)) {
            LOGGER.warn("Could not find gene {} for event {} in HMF driver gene panel. Skipping extraction!", gene, event);
        }
    }

    public static boolean checkGensInPanelForCuration(@NotNull String gene, @NotNull String event) {
        if (!GENES.contains(gene)) {
            LOGGER.warn("Could not find gene {} for event {} in HMF driver gene panel. Skipping extraction!", gene, event);
            return false;
        } else {
            return true;
        }
    }

}
