package com.hartwig.hmftools.orange.algo.sage;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.metrics.GeneDepth;
import com.hartwig.hmftools.common.metrics.GeneDepthFile;

import org.jetbrains.annotations.Nullable;

public final class GermlineMVLHFactory
{
    public static Map<String, Double> loadGermlineMVLHPerGene(
            final String sageGermlineGeneCoverageTsv,
            final Map<String,DriverGene> driverGenes) throws IOException
    {
        List<GeneDepth> geneDepths = GeneDepthFile.read(sageGermlineGeneCoverageTsv);
        return parseMVLHPerGene(geneDepths, driverGenes);
    }

    @VisibleForTesting
    static Map<String, Double> parseMVLHPerGene(final List<GeneDepth> geneDepths, final Map<String,DriverGene> driverGenes)
    {
        Map<String, Double> mvlhPerGene = Maps.newTreeMap();

        for(GeneDepth geneDepth : geneDepths)
        {
            DriverGene driverGene = driverGenes.get(geneDepth.Gene);

            if(driverGene != null && hasReliableMVLH(driverGene))
            {
                mvlhPerGene.put(geneDepth.Gene, geneDepth.MissedVariantLikelihood);
            }
        }
        return mvlhPerGene;
    }

    private static boolean hasReliableMVLH(final DriverGene driverGene)
    {
        // Note: fix this when SAGE germline starts running genome wide
        return driverGene.reportSomatic() || driverGene.reportGermline() || driverGene.reportPGX();
    }
}
