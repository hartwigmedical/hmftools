package com.hartwig.hmftools.peach;

import com.google.common.io.Resources;
import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

import org.jetbrains.annotations.NotNull;

public class TestUtils
{
    @NotNull
    public static HaplotypePanel loadTestHaplotypePanel(@NotNull String fileName)
    {
        return PanelLoader.loadHaplotypePanel(getTestResourcePath(fileName));
    }

    @NotNull
    public static String getTestResourcePath(@NotNull String fileName)
    {
        return Resources.getResource(fileName).getPath();
    }
}
