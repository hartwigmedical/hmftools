package com.hartwig.hmftools.serve.checkertool;

import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CheckExons {
    private static final Logger LOGGER = LogManager.getLogger(CheckExons.class);

    public CheckExons() {

    }

    public void run(@NotNull String event, @NotNull String geneSymbol) {
        LOGGER.info("Event check exons {} of gene {}", event, geneSymbol);
        Map<String, String> exonMap = Maps.newHashMap();
        exonMap.put(geneSymbol, event);

    }
}
