package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cli.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotFile;

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
    String SOMATIC_HOTSPOT = "somatic_hotspots";
    String GERMLINE_HOTSPOT = "germline_hotspots";

    static void addOptions(@NotNull Options options) {
        options.addOption(DRIVER_ENABLED, false, "Persist data to DB.");
        options.addOption(SOMATIC_HOTSPOT, true, "Path to somatic hotspot VCF");
        options.addOption(GERMLINE_HOTSPOT, true, "Path to germline hotspot VCF");
        DriverGenePanelConfig.addGenePanelOption(false, options);
    }

    boolean enabled();

    @NotNull
    DriverGenePanel genePanel();

    @NotNull
    ListMultimap<Chromosome, VariantHotspot> somaticHotspots();

    @NotNull
    ListMultimap<Chromosome, VariantHotspot> germlineHotspots();

    @NotNull
    static DriverCatalogConfig createConfig(@NotNull final CommandLine cmd, @NotNull RefGenomeData refGenomeData, @NotNull GermlineConfig germlineConfig)
            throws ParseException, IOException {
        boolean enabled = cmd.hasOption(DRIVER_ENABLED);
        String somaticHotspotVcf = cmd.getOptionValue(SOMATIC_HOTSPOT, Strings.EMPTY);
        String germlineHotspotVcf = cmd.getOptionValue(GERMLINE_HOTSPOT, Strings.EMPTY);
        final DriverGenePanel genePanel;

        if (enabled) {
            if (!DriverGenePanelConfig.isConfigured(cmd)) {
                throw new ParseException(
                        DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION + " is a mandatory argument when " + DRIVER_ENABLED + " enabled");
            }

            if (somaticHotspotVcf.isEmpty()) {
                throw new ParseException(SOMATIC_HOTSPOT + " is a mandatory argument when " + DRIVER_ENABLED + " enabled");
            }

            if (!new File(somaticHotspotVcf).exists()) {
                throw new IOException("Unable to open " + SOMATIC_HOTSPOT + " file " + somaticHotspotVcf);
            }

            final List<DriverGene> driverGenes = DriverGenePanelConfig.driverGenes(cmd);
            final RefGenomeVersion driverGeneRefGenomeVersion =
                    refGenomeData.isHg38() ? RefGenomeVersion.V38 : RefGenomeVersion.V37;
            genePanel = DriverGenePanelFactory.create(driverGeneRefGenomeVersion, driverGenes);

            if (germlineConfig.enabled()) {
                if (germlineHotspotVcf.isEmpty()) {
                    throw new ParseException(GERMLINE_HOTSPOT + " is a mandatory argument when " + DRIVER_ENABLED + " enabled");
                }

                if (!new File(germlineHotspotVcf).exists()) {
                    throw new IOException("Unable to open " + GERMLINE_HOTSPOT + " file " + germlineHotspotVcf);
                }
            }

        } else {
            genePanel = DriverGenePanelFactory.empty();
        }

        ListMultimap<Chromosome, VariantHotspot> somaticHotspots =
                somaticHotspotVcf.equals(Strings.EMPTY) ? ArrayListMultimap.create() : VariantHotspotFile.readFromVCF(somaticHotspotVcf);

        ListMultimap<Chromosome, VariantHotspot> germlineHotspots =
                germlineHotspotVcf.equals(Strings.EMPTY) ? ArrayListMultimap.create() : VariantHotspotFile.readFromVCF(germlineHotspotVcf);

        return ImmutableDriverCatalogConfig.builder()
                .enabled(enabled)
                .somaticHotspots(somaticHotspots)
                .germlineHotspots(germlineHotspots)
                .genePanel(genePanel)
                .build();
    }
}
