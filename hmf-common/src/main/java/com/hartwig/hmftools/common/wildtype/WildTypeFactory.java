package com.hartwig.hmftools.common.wildtype;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
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
            @NotNull List<DriverGene> driverGenes) {

        List<WildType> wildTypeGenes = Lists.newArrayList();
        boolean wildTypeSomaticVariant = false;
        boolean wildTypeGermlineVariant = false;
        boolean wildTypeSomaticGainLoss = false;
        boolean wildTypeFusions = false;
        boolean wildTypeHomozygousDisruption = false;

        List<Boolean> somaticVariantBoolean = Lists.newArrayList();
        List<Boolean> somaticGermlineBoolean = Lists.newArrayList();
        List<Boolean> CNVBoolean = Lists.newArrayList();
        List<Boolean> FusionBooleanEnd = Lists.newArrayList();
        List<Boolean> FusionBooleanStart = Lists.newArrayList();
        List<Boolean> HomozygousDisruptionBoolean = Lists.newArrayList();

        for (DriverGene driverGene : driverGenes) {

            for (ReportableVariant somaticVariant : reportableSomaticVariant) {
                LOGGER.info("driverGene: " + driverGene);
                isWildType(driverGene.gene(), somaticVariant.gene(), somaticVariantBoolean);
                wildTypeSomaticVariant = somaticVariantBoolean.contains(true);
                LOGGER.info("wildTypeSomaticVariant: " + wildTypeSomaticVariant);
            }

            for (ReportableVariant germlineVariant : reportableGermlineVariant) {
                isWildType(driverGene.gene(), germlineVariant.gene(), somaticGermlineBoolean);
                wildTypeGermlineVariant = somaticGermlineBoolean.contains(true);
                LOGGER.info("wildTypeGermlineVariant: " + wildTypeGermlineVariant);

            }

            for (GainLoss gainLoss : reportableSomaticGainsLosses) {
                isWildType(driverGene.gene(), gainLoss.gene(), CNVBoolean);
                wildTypeSomaticGainLoss = CNVBoolean.contains(true);
                LOGGER.info("wildTypeSomaticGainLoss: " + wildTypeSomaticGainLoss);

            }

            for (LinxFusion fusions : reportableFusions) {
                isWildType(driverGene.gene(), fusions.geneStart(), FusionBooleanStart);
                isWildType(driverGene.gene(), fusions.geneEnd(), FusionBooleanEnd);
                LOGGER.info("FusionBooleanStart: " + FusionBooleanStart);
                LOGGER.info("FusionBooleanEnd: " + FusionBooleanEnd);

                wildTypeFusions =
                        (!FusionBooleanStart.contains(true) && !FusionBooleanEnd.contains(true)) || (FusionBooleanStart.contains(true)
                                && FusionBooleanEnd.contains(true));
                LOGGER.info("wildTypeFusions: " + wildTypeFusions);

            }

            for (ReportableHomozygousDisruption homozygousDisruption : homozygousDisruptions) {

                isWildType(driverGene.gene(), homozygousDisruption.gene(), HomozygousDisruptionBoolean);
                wildTypeHomozygousDisruption = HomozygousDisruptionBoolean.contains(true);
                LOGGER.info("wildTypeHomozygousDisruption: " + wildTypeHomozygousDisruption);

            }

            if (!wildTypeSomaticVariant && !wildTypeGermlineVariant && !wildTypeSomaticGainLoss && !wildTypeFusions
                    && !wildTypeHomozygousDisruption) {
                wildTypeGenes.add(ImmutableWildType.builder().gene(driverGene.gene()).build());
            }
        }
        return wildTypeGenes;
    }

    public static void isWildType(@NotNull String driverGene, @NotNull String reportableGene, @NotNull List<Boolean> booleanList) {
        LOGGER.info("driverGene: " + driverGene );
        LOGGER.info("reportableGene: " + reportableGene);
        if (driverGene.equals(reportableGene)) {
            booleanList.add(true);
        } else {
            booleanList.add(false);
        }
    }
}