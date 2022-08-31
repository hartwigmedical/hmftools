package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public final class DriverGenePanelConfig
{
    public static final String DRIVER_GENE_PANEL_OPTION = "driver_gene_panel";
    public static final String DRIVER_GENE_PANEL_OPTION_DESC = "Path to driver gene panel";

    private DriverGenePanelConfig() {}

    public static void addGenePanelOption(boolean isRequired, final Options options)
    {
        Option genePanelOption = new Option(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);
        genePanelOption.setRequired(isRequired);
        options.addOption(genePanelOption);
    }

    public static boolean isConfigured(@NotNull final CommandLine cmd)
    {
        return cmd.hasOption(DRIVER_GENE_PANEL_OPTION);
    }

    @NotNull
    public static List<DriverGene> driverGenes(@NotNull final CommandLine cmd) throws IOException
    {
        return DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));
    }

    public static List<DriverGene> loadDriverGenes(final CommandLine cmd)
    {
        if(!isConfigured(cmd))
            return Lists.newArrayList();

        try
        {
            return DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));
        }
        catch(IOException e)
        {
            return Lists.newArrayList();
        }
    }
}
