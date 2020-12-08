package com.hartwig.hmftools.serve.util;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProgressTracker {

    private static final Logger LOGGER = LogManager.getLogger(ProgressTracker.class);
    private static final int INTERVAL = 10;

    @NotNull
    private final String label;
    private final int totalCount;

    private int counter = 0;

    public ProgressTracker(@NotNull final String label, final int totalCount) {
        this.label = label;
        this.totalCount = totalCount;
    }

    public void update() {
        ++counter;
        if (totalCount >= INTERVAL && counter % (totalCount / INTERVAL) == 0) {
            LOGGER.info(" Processed {} of {} {} entries", counter, totalCount, label);
        }
    }
}
