package com.hartwig.hmftools.orange.algo.wildtype;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.loader.GainLoss;
import com.hartwig.hmftools.orange.algo.linx.GeneDisruption;
import com.hartwig.hmftools.orange.algo.purple.ReportableVariant;

import org.jetbrains.annotations.NotNull;

public final class WildTypeFactory {

    private WildTypeFactory() {
    }

    @NotNull
    public static List<WildTypeGene> filterQCWildTypes(@NotNull Set<PurpleQCStatus> purpleQCStatus,
            @NotNull List<WildTypeGene> wildTypeGenes) {
        if (!purpleQCStatus.contains(PurpleQCStatus.FAIL_NO_TUMOR) && !purpleQCStatus.contains(PurpleQCStatus.WARN_LOW_PURITY)) {
            return wildTypeGenes;
        }
        return Lists.newArrayList();
    }

    @NotNull
    public static List<WildTypeGene> determineWildTypeGenes(@NotNull List<ReportableVariant> reportableGermlineVariants,
            @NotNull List<ReportableVariant> reportableSomaticVariants, @NotNull List<GainLoss> reportableSomaticGainsLosses,
            @NotNull List<LinxFusion> reportableFusions, @NotNull List<HomozygousDisruption> homozygousDisruptions,
            @NotNull List<GeneDisruption> reportableGeneDisruptions, @NotNull List<DriverGene> driverGenes) {
        List<WildTypeGene> wildTypeGenes = Lists.newArrayList();

        for (DriverGene driverGene : driverGenes) {
            boolean hasSomaticVariant = false;
            for (ReportableVariant somaticVariant : reportableSomaticVariants) {
                if (driverGene.gene().equals(somaticVariant.gene())) {
                    hasSomaticVariant = true;
                }
            }

            boolean hasGermlineVariant = false;
            for (ReportableVariant germlineVariant : reportableGermlineVariants) {
                if (driverGene.gene().equals(germlineVariant.gene())) {
                    hasGermlineVariant = true;
                }
            }

            boolean hasSomaticGainLoss = false;
            for (GainLoss gainLoss : reportableSomaticGainsLosses) {
                if (driverGene.gene().equals(gainLoss.gene())) {
                    hasSomaticGainLoss = true;
                }
            }

            boolean hasFusion = false;
            for (LinxFusion fusion : reportableFusions) {
                if (driverGene.gene().equals(fusion.geneStart()) || driverGene.gene().equals(fusion.geneEnd())) {
                    hasFusion = true;
                }
            }

            boolean hasHomozygousDisruption = false;
            for (HomozygousDisruption homozygousDisruption : homozygousDisruptions) {
                if (driverGene.gene().equals(homozygousDisruption.gene())) {
                    hasFusion = true;
                }
            }

            boolean hasGeneDisruption = false;
            for (GeneDisruption disruption : reportableGeneDisruptions) {
                if (driverGene.gene().equals(disruption.gene())) {
                    hasGeneDisruption = true;
                }
            }

            if (!hasSomaticVariant && !hasGermlineVariant && !hasSomaticGainLoss && !hasFusion && !hasHomozygousDisruption
                    && !hasGeneDisruption) {
                wildTypeGenes.add(ImmutableWildTypeGene.builder().gene(driverGene.gene()).build());
            }
        }

        return wildTypeGenes;
    }
}