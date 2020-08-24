package com.hartwig.hmftools.common.cli;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class DriverGenePanelConfig {

    public static final String DRIVER_GENE_PANEL_OPTION = "driver_gene_panel";

    public static void addGenePanelOption(boolean isRequired, final Options options) {
        Option genePanelOption = new Option(DRIVER_GENE_PANEL_OPTION, true, "Path to driver gene panel");
        genePanelOption.setRequired(isRequired);
        options.addOption(genePanelOption);
    }

    public static boolean isConfigured(@NotNull final CommandLine cmd) {
        return cmd.hasOption(DRIVER_GENE_PANEL_OPTION);
    }

    @NotNull
    public static DriverGenePanel driverGenePanel(@NotNull final CommandLine cmd) throws IOException {
        List<DriverGene> driverGenes = DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));
        return new DriverGenePanelFactory().create(driverGenes);
    }

}
