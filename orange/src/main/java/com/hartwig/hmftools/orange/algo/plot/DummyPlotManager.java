package com.hartwig.hmftools.orange.algo.plot;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class DummyPlotManager implements PlotManager {

    private static final Logger LOGGER = LogManager.getLogger(DummyPlotManager.class);

    @Override
    public void createPlotDirectory() {
        LOGGER.debug("Creating dummy plot directory");
    }

    @Nullable
    @Override
    public String processPlotFile(@Nullable final String sourcePlotPath) {
        LOGGER.debug("Dummy-processing plot file '{}'", sourcePlotPath);
        return sourcePlotPath;
    }
}
