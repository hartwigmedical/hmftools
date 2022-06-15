package com.hartwig.hmftools.common.wildtype;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class WildTypeFactory {
    private static final Logger LOGGER = LogManager.getLogger(WildTypeFactory.class);

    public WildTypeFactory() {
    }

    public static List<WildType> determineWildTypeGenes(@NotNull List<ReportableVariant> reportableGermlineVariant,
            @NotNull List<ReportableVariant> reportableSomaticVariant, @NotNull List<GainLoss> reportableSomaticGainsLosses,
            @NotNull List<LinxFusion> reportableFusions, @NotNull List<ReportableHomozygousDisruption> homozygousDisruptions,
            @NotNull List<ReportableGeneDisruption> geneDisruptions, @NotNull List<DriverGene> driverGenes) {

        List<WildType> wildTypeGenes = Lists.newArrayList();

        for (DriverGene driverGene : driverGenes) {
            LOGGER.info("driverGene: " + driverGene);

            boolean hasSomaticVariant = false;
            boolean hasGermlineVariant = false;
            boolean hasSomaticGainLoss = false;
            boolean hasFusion = false;
            boolean hasHomozygousDisruption = false;
            boolean hasGeneDisruption = false;

            for (ReportableVariant somaticVariant : reportableSomaticVariant) {
                if (driverGene.gene().equals(somaticVariant.gene())) {
                    hasSomaticVariant = true;
                }
            }
            LOGGER.info("hasSomaticVariant: " + hasSomaticVariant);

            for (ReportableVariant germlineVariant : reportableGermlineVariant) {
                if (driverGene.gene().equals(germlineVariant.gene())) {
                    hasGermlineVariant = true;
                }
            }
            LOGGER.info("hasGermlineVariant: " + hasGermlineVariant);

            for (GainLoss gainLoss : reportableSomaticGainsLosses) {
                if (driverGene.gene().equals(gainLoss.gene())) {
                    hasSomaticGainLoss = true;
                }
            }
            LOGGER.info("hasSomaticGainLoss: " + hasSomaticGainLoss);

            for (LinxFusion fusion : reportableFusions) {
                if (driverGene.gene().equals(fusion.geneStart())) {
                    hasFusion = true;
                }

                if (driverGene.gene().equals(fusion.geneEnd())) {
                    hasFusion = true;
                }
            }
            LOGGER.info("hasFusion: " + hasFusion);

            for (ReportableHomozygousDisruption homozygousDisruption : homozygousDisruptions) {
                if (driverGene.gene().equals(homozygousDisruption.gene())) {
                    hasFusion = true;
                }
            }
            LOGGER.info("hasHomozygousDisruption: " + hasHomozygousDisruption);

            for (ReportableGeneDisruption geneDisruption : geneDisruptions) {
                LOGGER.info(geneDisruption);
                LOGGER.info(geneDisruption.gene());
                LOGGER.info(geneDisruption.gene());
                if (driverGene.gene().equals(geneDisruption.gene())) {
                    hasGeneDisruption = true;
                }
            }
            LOGGER.info("hasGeneDisruption: " + hasGeneDisruption);


            if (!hasSomaticVariant && !hasGermlineVariant && !hasSomaticGainLoss && !hasFusion && !hasHomozygousDisruption
                    && !hasGeneDisruption) {
                wildTypeGenes.add(ImmutableWildType.builder().gene(driverGene.gene()).build());
            }
        }
        return wildTypeGenes;
    }
}