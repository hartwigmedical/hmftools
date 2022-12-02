package com.hartwig.hmftools.orange.algo.plot;

import java.io.IOException;

import org.jetbrains.annotations.Nullable;

public interface PlotManager {

    void createPlotDirectory() throws IOException;

    @Nullable
    String processPlotFile(@Nullable String sourcePlotPath) throws IOException;
}
