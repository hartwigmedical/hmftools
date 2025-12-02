package com.hartwig.hmftools.orange.algo.wildtype;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.datamodel.finding.Fusion;
import com.hartwig.hmftools.datamodel.finding.GainDeletion;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.wildtype.ImmutableWildTypeGene;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;

import org.jetbrains.annotations.Nullable;

public final class WildTypeAlgo
{
    public static boolean wildTypeCallingAllowed(final Set<PurpleQCStatus> purpleQCStatus)
    {
        return !purpleQCStatus.contains(PurpleQCStatus.FAIL_NO_TUMOR) && !purpleQCStatus.contains(PurpleQCStatus.WARN_LOW_PURITY);
    }

    public static List<WildTypeGene> determineWildTypeGenes(
            final Map<String,DriverGene> driverGenes,
            final List<PurpleVariant> reportableSomaticVariants, @Nullable List<PurpleVariant> reportableGermlineVariants,
            final List<GainDeletion> reportableSomaticGainsDels, final List<Fusion> reportableFusions,
            final List<LinxHomozygousDisruption> homozygousDisruptions, final List<LinxBreakend> reportableBreakends)
    {
        List<WildTypeGene> wildTypeGenes = Lists.newArrayList();

        for(DriverGene driverGene : driverGenes.values())
        {
            boolean hasSomaticVariant = false;
            for(PurpleVariant somaticVariant : reportableSomaticVariants)
            {
                if(driverGene.gene().equals(somaticVariant.gene()))
                {
                    hasSomaticVariant = true;
                }
            }

            boolean hasGermlineVariant = false;
            if(reportableGermlineVariants != null)
            {
                for(PurpleVariant germlineVariant : reportableGermlineVariants)
                {
                    if(driverGene.gene().equals(germlineVariant.gene()))
                    {
                        hasGermlineVariant = true;
                    }
                }
            }

            boolean hasSomaticgainDel = false;
            for(GainDeletion gainDel : reportableSomaticGainsDels)
            {
                if(driverGene.gene().equals(gainDel.gene()))
                {
                    hasSomaticgainDel = true;
                }
            }

            boolean hasFusion = false;
            for(Fusion fusion : reportableFusions)
            {
                if(driverGene.gene().equals(fusion.geneStart()) || driverGene.gene().equals(fusion.geneEnd()))
                {
                    hasFusion = true;
                }
            }

            boolean hasHomozygousDisruption = false;
            for(LinxHomozygousDisruption homozygousDisruption : homozygousDisruptions)
            {
                if(driverGene.gene().equals(homozygousDisruption.gene()))
                {
                    hasHomozygousDisruption = true;
                }
            }

            boolean hasBreakend = false;
            for(LinxBreakend breakend : reportableBreakends)
            {
                if(driverGene.gene().equals(breakend.gene()))
                {
                    hasBreakend = true;
                }
            }

            if(!hasSomaticVariant && !hasGermlineVariant && !hasSomaticgainDel && !hasFusion && !hasHomozygousDisruption && !hasBreakend)
            {
                wildTypeGenes.add(ImmutableWildTypeGene.builder().gene(driverGene.gene()).build());
            }
        }

        return wildTypeGenes;
    }
}