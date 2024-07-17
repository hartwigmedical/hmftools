package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.FusionPhasedType;
import com.hartwig.hmftools.common.linx.LinxFusion;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class DnaFusionSelector
{
    @NotNull
    public static List<LinxFusion> selectInterestingUnreportedFusions(@NotNull List<LinxFusion> allFusions,
            @NotNull List<DriverGene> driverGenes)
    {
        List<LinxFusion> filtered = Lists.newArrayList();
        for(LinxFusion fusion : allFusions)
        {
            if(!fusion.reported())
            {
                boolean hasReportedType = !fusion.reportedType().equals(KnownFusionType.NONE.toString());
                boolean isFusionOfOncogene = isInframeFusionWithOncogene(fusion, driverGenes);
                if(hasReportedType || isFusionOfOncogene)
                {
                    filtered.add(fusion);
                }
            }
        }
        return filtered;
    }

    @NotNull
    public static List<LinxFusion> selectAdditionalViableSomaticFusions(@NotNull List<LinxFusion> somaticFusions,
            @NotNull List<LinxFusion> additionalSuspectSomaticFusions)
    {
        return somaticFusions.stream()
                .filter(fusion -> fusion.phased() == FusionPhasedType.INFRAME && !fusion.chainTerminated() && !fusion.reported()
                        && !additionalSuspectSomaticFusions.contains(fusion))
                .collect(Collectors.toList());
    }

    private static boolean isInframeFusionWithOncogene(@NotNull LinxFusion fusion, @NotNull List<DriverGene> driverGenes)
    {
        if(fusion.phased() != FusionPhasedType.INFRAME)
        {
            return false;
        }

        return isOncoDriverGene(driverGenes, fusion.geneStart()) || isOncoDriverGene(driverGenes, fusion.geneEnd());
    }

    private static boolean isOncoDriverGene(@NotNull List<DriverGene> driverGenes, @NotNull String gene)
    {
        DriverGene driver = findDriverGene(driverGenes, gene);
        return driver != null && driver.likelihoodType() == DriverCategory.ONCO;
    }

    @Nullable
    private static DriverGene findDriverGene(@NotNull List<DriverGene> driverGenes, @NotNull String geneToFind)
    {
        for(DriverGene driverGene : driverGenes)
        {
            if(driverGene.gene().equals(geneToFind))
            {
                return driverGene;
            }
        }
        return null;
    }
}
