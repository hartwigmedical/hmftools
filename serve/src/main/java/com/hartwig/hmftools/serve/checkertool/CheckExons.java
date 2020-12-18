package com.hartwig.hmftools.serve.checkertool;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CheckExons {
    private static final Logger LOGGER = LogManager.getLogger(CheckExons.class);

    public CheckExons() {

    }

    public void run(@NotNull String event) {
        LOGGER.info("Event check exons: ", event);

    }
}
