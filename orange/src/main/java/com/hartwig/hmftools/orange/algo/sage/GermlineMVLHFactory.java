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
            final List<DriverGene> driverGenes) throws IOException
    {
        List<GeneDepth> geneDepths = GeneDepthFile.read(sageGermlineGeneCoverageTsv);
        return parseMVLHPerGene(geneDepths, driverGenes);
    }

    @VisibleForTesting
    static Map<String, Double> parseMVLHPerGene(final List<GeneDepth> geneDepths, final List<DriverGene> driverGenes)
    {
        Map<String, Double> mvlhPerGene = Maps.newTreeMap();

        for(GeneDepth geneDepth : geneDepths)
        {
            DriverGene matchingDriverGene = findMatchingDriverGene(geneDepth.Gene, driverGenes);
            if(matchingDriverGene != null && hasReliableMVLH(matchingDriverGene))
            {
                mvlhPerGene.put(geneDepth.Gene, geneDepth.MissedVariantLikelihood);
            }
        }
        return mvlhPerGene;
    }

    @Nullable
    private static DriverGene findMatchingDriverGene(final String geneName, final List<DriverGene> driverGenes)
    {
        List<DriverGene> matchingDriverGenes = driverGenes.stream().filter(d -> d.gene().equals(geneName)).collect(Collectors.toList());
        if(matchingDriverGenes.size() == 1)
        {
            return matchingDriverGenes.get(0);
        }
        else if(matchingDriverGenes.isEmpty())
        {
            return null;
        }
        else
        {
            throw new IllegalStateException("More than one driver gene found with name " + geneName);
        }
    }

    private static boolean hasReliableMVLH(final DriverGene driverGene)
    {
        // Note: fix this when SAGE germline starts running genome wide
        return driverGene.reportSomatic() || driverGene.reportGermline() || driverGene.reportPGX();
    }
}
