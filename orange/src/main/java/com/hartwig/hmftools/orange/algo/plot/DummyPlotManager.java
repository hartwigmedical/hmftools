package com.hartwig.hmftools.orange.algo.plot;

import org.jetbrains.annotations.Nullable;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

public class DummyPlotManager implements PlotManager
{
    @Override
    public void createPlotDirectory()
    {
        LOGGER.debug("Creating dummy plot directory");
    }

    @Nullable
    @Override
    public String processPlotFile(@Nullable final String sourcePlotPath)
    {
        LOGGER.debug("Dummy-processing plot file '{}'", sourcePlotPath);
        return sourcePlotPath;
    }
}
