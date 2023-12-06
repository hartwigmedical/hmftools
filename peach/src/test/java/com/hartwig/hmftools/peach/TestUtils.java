package com.hartwig.hmftools.peach;

import com.google.common.io.Resources;
import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

public class TestUtils
{
    public static HaplotypePanel loadTestHaplotypePanel(String fileName)
    {
        return PanelLoader.loadHaplotypePanel(getTestResourcePath(fileName));
    }

    public static String getTestResourcePath(String fileName)
    {
        return Resources.getResource(fileName).getPath();
    }
}
