package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface DriverCatalogConfig {

    String DRIVER_ENABLED = "driver_catalog";
    String HOTSPOT = "hotspots";
    String DRIVER_GENE_PANEL = "gene_panel";

    static void addOptions(@NotNull Options options) {
        options.addOption(DRIVER_ENABLED, false, "Persist data to DB.");
        options.addOption(HOTSPOT, true, "Database user name.");
        options.addOption(DRIVER_GENE_PANEL, true, "Driver gene panel.");
    }

    boolean enabled();

    @NotNull
    String hotspots();

    @NotNull
    DriverGenePanel genePanel();

    @NotNull
    static DriverCatalogConfig createConfig(@NotNull final CommandLine cmd) throws ParseException, IOException {
        boolean enabled = cmd.hasOption(DRIVER_ENABLED);
        String hotspots = cmd.getOptionValue(HOTSPOT, Strings.EMPTY);
        DriverGenePanel genePanel;
        if (cmd.hasOption(DRIVER_GENE_PANEL)) {
            List<DriverGene> driverGenes = DriverGeneFile.read(cmd.getOptionValue(DRIVER_GENE_PANEL));
            genePanel = new DriverGenePanelFactory().create(driverGenes);
        } else {
            genePanel = new DriverGenePanelFactory().create();
        }

        if (enabled) {
            if (hotspots.isEmpty()) {
                throw new ParseException(HOTSPOT + " is a mandatory argument when " + DRIVER_ENABLED + " enabled");
            }

            if (!new File(hotspots).exists()) {
                throw new IOException("Unable to open " + HOTSPOT + " file " + hotspots);
            }

        }
        return ImmutableDriverCatalogConfig.builder().enabled(enabled).hotspots(hotspots).genePanel(genePanel).build();
    }
}
