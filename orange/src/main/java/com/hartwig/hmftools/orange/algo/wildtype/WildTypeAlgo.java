package com.hartwig.hmftools.orange.algo.wildtype;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class WildTypeAlgo {

    private WildTypeAlgo() {
    }

    public static boolean wildTypeCallingAllowed(@NotNull Set<PurpleQCStatus> purpleQCStatus) {
        return !purpleQCStatus.contains(PurpleQCStatus.FAIL_NO_TUMOR) && !purpleQCStatus.contains(PurpleQCStatus.WARN_LOW_PURITY);
    }

    @NotNull
    public static List<WildTypeGene> determineWildTypeGenes(@NotNull List<DriverGene> driverGenes,
            @NotNull List<PurpleVariant> reportableSomaticVariants, @Nullable List<PurpleVariant> reportableGermlineVariants,
            @NotNull List<PurpleGainLoss> reportableSomaticGainsLosses, @NotNull List<LinxFusion> reportableFusions,
            @NotNull List<HomozygousDisruption> homozygousDisruptions, @NotNull List<LinxBreakend> reportableBreakends) {
        List<WildTypeGene> wildTypeGenes = Lists.newArrayList();

        for (DriverGene driverGene : driverGenes) {
            boolean hasSomaticVariant = false;
            for (PurpleVariant somaticVariant : reportableSomaticVariants) {
                if (driverGene.gene().equals(somaticVariant.gene())) {
                    hasSomaticVariant = true;
                }
            }

            boolean hasGermlineVariant = false;
            if (reportableGermlineVariants != null) {
                for (PurpleVariant germlineVariant : reportableGermlineVariants) {
                    if (driverGene.gene().equals(germlineVariant.gene())) {
                        hasGermlineVariant = true;
                    }
                }
            }

            boolean hasSomaticGainLoss = false;
            for (PurpleGainLoss gainLoss : reportableSomaticGainsLosses) {
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

            boolean hasBreakend = false;
            for (LinxBreakend breakend : reportableBreakends) {
                if (driverGene.gene().equals(breakend.gene())) {
                    hasBreakend = true;
                }
            }

            if (!hasSomaticVariant && !hasGermlineVariant && !hasSomaticGainLoss && !hasFusion && !hasHomozygousDisruption
                    && !hasBreakend) {
                wildTypeGenes.add(ImmutableWildTypeGene.builder().gene(driverGene.gene()).build());
            }
        }

        return wildTypeGenes;
    }
}