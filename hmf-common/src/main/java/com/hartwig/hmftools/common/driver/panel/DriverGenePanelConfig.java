package com.hartwig.hmftools.common.driver.panel;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public final class DriverGenePanelConfig
{
    public static final String DRIVER_GENE_PANEL = "driver_gene_panel";
    public static final String DRIVER_GENE_PANEL_DESC = "Path to driver gene panel";

    private DriverGenePanelConfig() {}

    public static boolean isConfigured(final ConfigBuilder configBuilder)
    {
        return configBuilder.hasValue(DRIVER_GENE_PANEL);
    }

    public static List<DriverGene> loadDriverGenes(final ConfigBuilder configBuilder)
    {
        if(!isConfigured(configBuilder))
            return Lists.newArrayList();

        try
        {
            return DriverGeneFile.read(configBuilder.getValue(DRIVER_GENE_PANEL));
        }
        catch(IOException e)
        {
            return Lists.newArrayList();
        }
    }

    public static void addGenePanelOption(final ConfigBuilder configBuilder, boolean isRequired)
    {
        configBuilder.addPath(DRIVER_GENE_PANEL, isRequired, DRIVER_GENE_PANEL_DESC);
    }

    public static List<DriverGene> driverGenes(final ConfigBuilder configBuilder) throws IOException
    {
        return DriverGeneFile.read(configBuilder.getValue(DRIVER_GENE_PANEL));
    }
}
